#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-05-13 21:46:25
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-05-19 13:00:34

import scanpy as sc
import numpy as np
import pandas as pd
import scipy.stats
import scrublet as scr
from ...utilities._misc import *

def run_scanpy_qc(self, max_genes=2500, min_genes=200, mito_cutoff=0.05, pval_cutoff=0.1, min_counts=None, max_counts=None):
    """
    Parameters
    ----------
    adata : AnnData
        The (annotated) data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    max_genes : int (default: 2500)
        Maximum number of genes expressed required for a cell to pass filtering.
    min_genes : int (default: 200)
        Minimum number of genes expressed  required for a cell to pass filtering.
    mito_cutoff : float (default: 0.05)
        Maximum percentage mitochondrial content allowed for a cell to pass filtering.
    pval_cutoff : float (default: 0.05)
        Maximum Benjamini-Hochberg corrected p value from doublet detection protocol allowed for a cell to pass filtering.
    min_counts : int, None (default: None)
        Minimum number of counts required for a cell to pass filtering.
    max_counts : int, None (default: None)
        Maximum number of counts required for a cell to pass filtering.
    Returns
    -------
    AnnData
        The (annotated) data matrix of shape n_obs × n_vars where obs now contain filtering information. Rows correspond to cells and columns to genes.

    """
    _adata = self.copy()
    # run scrublet
    scrub = scr.Scrublet(_adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    _adata.obs['scrublet_score'] = doublet_scores
    # overcluster prep. run basic scanpy pipeline
    sc.pp.filter_cells(_adata, min_genes = 0)
    mito_genes = _adata.var_names.str.startswith('MT-')
    _adata.obs['percent_mito'] = np.sum(_adata[:, mito_genes].X, axis = 1) / np.sum(_adata.X, axis = 1)
    _adata.obs['n_counts'] = _adata.X.sum(axis = 1)
    sc.pp.normalize_total(_adata)
    sc.pp.log1p(_adata)
    sc.pp.highly_variable_genes(_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    _adata = _adata[:, _adata.var['highly_variable']]
    sc.pp.scale(_adata, max_value=10)
    sc.tl.pca(_adata, svd_solver='arpack')
    sc.pp.neighbors(_adata, n_neighbors=10, n_pcs=50)
    # overclustering proper - do basic clustering first, then cluster each cluster
    sc.tl.leiden(_adata)
    for clus in list(np.unique(_adata.obs['leiden']))[0]:
        sc.tl.leiden(_adata, restrict_to=('leiden',[clus]), key_added = 'leiden_R')
    for clus in list(np.unique(_adata.obs['leiden']))[1:]: # weird how the new anndata/scanpy is forcing this
        sc.tl.leiden(_adata, restrict_to=('leiden_R',[clus]), key_added = 'leiden_R')
    # compute the cluster scores - the median of Scrublet scores per overclustered cluster
    for clus in np.unique(_adata.obs['leiden_R']):
        _adata.obs.loc[_adata.obs['leiden_R']==clus, 'scrublet_cluster_score'] = \
            np.median(_adata.obs.loc[_adata.obs['leiden_R']==clus, 'scrublet_score'])
    # now compute doublet p-values. figure out the median and mad (from above-median values) for the distribution
    med = np.median(_adata.obs['scrublet_cluster_score'])
    mask = _adata.obs['scrublet_cluster_score']>med
    mad = np.median(_adata.obs['scrublet_cluster_score'][mask]-med)
    # let's do a one-sided test. the Bertie write-up does not address this but it makes sense
    pvals = 1-scipy.stats.norm.cdf(_adata.obs['scrublet_cluster_score'], loc=med, scale=1.4826*mad)
    _adata.obs['bh_pval'] = bh(pvals)
    # threshold the p-values to get doublet calls.
    _adata.obs['is_doublet'] = _adata.obs['bh_pval'] < pval_cutoff
    _adata.obs['is_doublet'] = _adata.obs['is_doublet'].astype('category')
    _adata.obs['filter_rna'] = (pd.Series([min_genes < n > max_genes for n in _adata.obs['n_genes']], index = _adata.obs.index)) | \
        (_adata.obs['percent_mito'] >= mito_cutoff) | \
            (_adata.obs['is_doublet'] == True)

    # removing columns that probably don't need anymore
    _adata.obs = _adata.obs.drop(['leiden', 'leiden_R', 'scrublet_cluster_score'], axis = 1)
    self.obs = _adata.obs.copy()
