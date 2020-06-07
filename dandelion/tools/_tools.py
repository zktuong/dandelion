#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-05-13 23:22:18
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-06-07 23:10:49

import os
import sys
import scanpy as sc
import pandas as pd
from pandas import DataFrame
import numpy as np
from tqdm import tqdm
from ..utilities._misc import *
from collections import defaultdict
from itertools import groupby
from scipy.spatial.distance import pdist, squareform
from distance import hamming
import re
import math
import scipy
from scipy.sparse.csgraph import minimum_spanning_tree
import Levenshtein
import networkx as nx
import igraph
from time import sleep
import copy
import functools
try:
    from scanpy import logging as logg
except ImportError:
    pass
from rpy2.robjects.packages import importr, data
from rpy2.rinterface import NULL
from rpy2.robjects import pandas2ri
import warnings
from subprocess import run
import multiprocessing
from changeo.Gene import getGene
from rpy2.robjects.packages import importr, data
from rpy2.rinterface import NULL
from rpy2.robjects import pandas2ri
import rpy2.robjects
from plotnine import ggplot, geom_point, options, annotate, aes, xlab, ylab, facet_grid, theme_bw, geom_histogram, geom_vline, facet_grid, theme


def filter_bcr(data, adata, filter_bcr=True, filter_rna=True, rescue_igh=True, umi_foldchange_cutoff=5, filter_lightchains=True, filter_missing=True, outdir=None, outFilePrefix=None, filtered=False):
    """
    Filters doublets and poor quality cells and corresponding contigs based on provided V(D)J `DataFrame` and `AnnData` objects. Depends on a `AnnData`.obs slot populated with 'filter_rna' column.
    Cells with multiple IGH contigs are filtered unless rescue_igh is True, where by the umi counts for each IGH contig will then be compared. The contig with the highest umi that is > umi_foldchange_cutoff (default is empirically set at 5) from the lowest will be retained.
    If there's multiple contigs that survive the 'rescue', then all contigs will be filtered. The default behaviour is to also filter cells with multiple lightchains but this may sometimes be a true biological occurrence; toggling filter_lightchains to False will rescue the mutltiplet light chains.
    Lastly, contigs with no corresponding cell barcode in the AnnData object is filtered if filter_missing is True. However, this may be useful to toggle to False if more contigs are preferred to be kept or for integrating with bulk reperotire seq data.

    Parameters
    ----------
    data : DataDrame, str
        V(D)J airr/changeo data to filter. Can be pandas `DataFrame` object or file path as string.
    adata : AnnData
        AnnData object to filter.
    filter_bcr : bool
        If True, V(D)J `DataFrame` object returned will be filtered. Default is True.
    filter_rna : bool
        If True, `AnnData` object returned will be filtered. Default is True.
    rescue_igh : bool
        If True, rescues IGH contigs with highest umi counts with a requirement that it passes the `umi_foldchange_cutoff` option. Default is True.
    umi_foldchange_cutoff : int
        related tominimum fold change required to rescue heavy chain contigs/barcode otherwise they will be marked as doublets. Default is empirically set at 5.
    filter_lightchains : bool
        cells with multiple light chains will be marked to filter. Default is True.
    filter_missing : bool
        cells in V(D)J data not found in `AnnData` object will be marked to filter. Default is True. This may be useful for toggling to False if integrating with bulk data.
    outdir : str, optional
        If specified, out file will be in this location
    outFilePrefix : str, optional
        If specified, the out file name will have this prefix
    filtered : bool
        If True, will create filenames with 'filtered_contig' as prefix. if False, will create filenames with 'all_contig' as prefix. ignored if outFilePrefix is specified.
    Returns
    -------
        V(D)J `DataFrame` object in airr/changeo format and `AnnData` object.
    """
    dat = load_data(data)
    h = Tree()
    h_umi = Tree()
    l = Tree()
    poor_qual = []
    h_doublet = []
    l_doublet = []
    drop_contig = []
    
    locus_dict = dict(zip(dat['sequence_id'],dat['locus']))
    barcode = list(set(dat['cell_id']))

    if 'filter_rna' not in adata.obs:
        raise TypeError("AnnData obs does not contain 'filter_rna' column. Please run `pp.recipe_scanpy_qc` before continuing.")

    bcr_check = Tree()
    for c in adata.obs_names:
        if c in barcode:
            bcr_check[c] = True
        else:
            bcr_check[c] = False
    adata.obs['has_bcr'] = pd.Series(dict(bcr_check))
    adata.obs['has_bcr'] = adata.obs['has_bcr'].astype('category')

    if 'v_call_genotyped' in dat.columns:
        v_dict = dict(zip(dat['sequence_id'], dat['v_call_genotyped']))
    else:
        v_dict = dict(zip(dat['sequence_id'], dat['v_call']))
    j_dict = dict(zip(dat['sequence_id'], dat['j_call']))

    for b in tqdm(barcode, desc = 'Marking barcodes with poor quality BCRs and BCR doublets'):
        hc_id = list(dat[(dat['cell_id'].isin([b])) & (dat['locus'] == 'IGH')]['sequence_id'])
        lc_id = list(dat[(dat['cell_id'].isin([b])) & (dat['locus'].isin(['IGK', 'IGL']))]['sequence_id'])
        hc_umi = [int(x) for x in dat[(dat['cell_id'].isin([b])) & (dat['locus'] == 'IGH')]['umi_count']]

        h[b] = hc_id
        h_umi[b] = hc_umi
        l[b] = lc_id
        # marking doublets defined by heavy chains
        if len(h[b]) > 1:
            if rescue_igh:
                highest_umi = max(h_umi[b])
                lowest_umi = min(h_umi[b])
                highest_umi_idx = [i for i, j in enumerate(h_umi[b]) if j == highest_umi]
                if len(highest_umi_idx) > 1:
                    h_doublet.append(b)
                if highest_umi/lowest_umi < umi_foldchange_cutoff:
                    h_doublet.append(b)
                if len(highest_umi_idx) == 1 and highest_umi/lowest_umi >= umi_foldchange_cutoff:
                    drop_contig.append(h[b][~highest_umi_idx[0]])
            else:
                h_doublet.append(b)
        # marking doublets defined by light chains
        if (len(h[b]) == 1) & (len(l[b]) > 1):
            l_doublet.append(b)
        # marking poor bcr quality, defined as those with only light chains, those
        # that were have conflicting assignment of locus and heavy/light V/J calls,
        # and also those that are missing either v or j calls
        if len(h[b]) < 1:
            poor_qual.append(b)
        if len(hc_id) > 0:
            v = v_dict[hc_id[0]]
            if 'IGH' not in v:
                poor_qual.append(b)
            j = j_dict[hc_id[0]]
            if 'IGH' not in j:
                poor_qual.append(b)
        if len(lc_id) > 0:
            v = v_dict[lc_id[0]]
            if 'IGH' in v:
                poor_qual.append(b)
            j = j_dict[lc_id[0]]
            if 'IGH' in j:
                poor_qual.append(b)

    poorqual = Tree()
    hdoublet = Tree()
    ldoublet = Tree()
    for c in tqdm(adata.obs_names, desc = 'Annotating in anndata obs slot '):
        if c in poor_qual:
            poorqual[c] = True
        else:
            poorqual[c] = False

        if c in h_doublet:
            hdoublet[c] = True
        else:
            hdoublet[c] = False

        if c in l_doublet:
            ldoublet[c] = True
        else:
            ldoublet[c] = False

    adata.obs['filter_bcr_quality'] = pd.Series(dict(poorqual))
    adata.obs['filter_bcr_quality'] = adata.obs['filter_bcr_quality'].astype('category')
    adata.obs['filter_bcr_heavy'] = pd.Series(dict(hdoublet))
    adata.obs['filter_bcr_heavy'] = adata.obs['filter_bcr_heavy'].astype('category')
    adata.obs['filter_bcr_light'] = pd.Series(dict(ldoublet))
    adata.obs['filter_bcr_light'] = adata.obs['filter_bcr_light'].astype('category')

    filter_ids = []
    if filter_bcr:
        if not filter_lightchains:
            filter_ids = list(set(h_doublet + poor_qual))
        else:
            filter_ids = list(set(h_doublet + l_doublet + poor_qual))

        if filter_rna:
            filter_ids = filter_ids + list(adata[adata.obs['filter_rna'] == True].obs_names)
            filter_ids = list(set(filter_ids))

        if filter_missing:
            for c in dat['cell_id']:
                if c not in adata.obs_names:
                    filter_ids.append(c)

        _dat = dat[~(dat['cell_id'].isin(filter_ids))]

        if rescue_igh:
            _dat = _dat[~(_dat['sequence_id'].isin(drop_contig))]

        if os.path.isfile(str(data)):
            _dat.to_csv("{}/{}_filtered.tsv".format(os.path.dirname(data), os.path.basename(data).split('.tsv')[0]), sep = '\t', index = None)
        else:
            if filtered:
                outFile_prefix = 'filtered_contig'
            else:
                outFile_prefix = 'all_contig'

            if (outdir is None) & (outFilePrefix is not None):
                _dat.to_csv("{}/{}_filtered.tsv".format('dandelion/data', str(outFilePrefix)), sep = '\t', index = None)
            elif (outdir is not None) & (outFilePrefix is None):
                _dat.to_csv("{}/{}_filtered.tsv".format(str(outdir), outFile_prefix), sep = '\t', index = None)
            elif (outdir is None) & (outFilePrefix is None):
                _dat.to_csv("{}/{}_filtered.tsv".format('dandelion/data', outFile_prefix), sep = '\t', index = None)
            elif (outdir is not None) & (outFilePrefix is None):
                _dat.to_csv("{}/{}_filtered.tsv".format(str(outdir), outFile_prefix), sep = '\t', index = None)

    if filter_rna:
        _adata = adata[~(adata.obs_names.isin(filter_ids))] # not saving the scanpy object because there's no need to at the moment
    else:
        _adata = adata.copy()

    return(_dat, _adata)

def find_clones(self, identity=0.85, clustering_by = None, by_alleles = None, write_out = False, outdir=None, outFilePrefix=None):
    """
    Find clones based on heavy chain and light chain CDR3 junction hamming distance.

    Parameters
    ----------
    self : Dandelion, DataFrame, str
        BCR data to find clone. Can be Dandelion object, pandas DataFrame or file path as string.
    identity : float
        Junction similarity parameter. Default 0.85    
    clustering_by : str, optional
        modes for clustering: 'nt' or 'aa'. None defaults to 'aa'.
    by_alleles : bool, optional
        Whether or not to collapse alleles to genes. None defaults to True.
    write_out : bool
        If True, will write out airr/changeo file with clone_id column (default is False). file path and file name is determined by outdir and outFilePrefix options.
    outdir : str, optional
        If specified, outfile will be in this location. None defaults to 'dandelion/data'.
    outFilePrefix : str, optional
        If specified, the outfile name will have this prefix. None defaults to 'dandelion_find_clones'
    Returns
    -------
        `Dandelion` object with clone_id annotated in `.data` slot and `.metadata` initialized.
    """
    start = logg.info('Finding clones')
    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    else:
        dat = load_data(self)
    dat_heavy = dat[dat['locus'] == 'IGH']
    pd.set_option('mode.chained_assignment', None)

    # retrieve the J genes and J genes
    if by_alleles is None or ~by_alleles:
        if 'v_call_genotyped' in dat_heavy.columns:
            V = [re.sub('[*][0-9][0-9]', '', v) for v in dat_heavy['v_call_genotyped']]
        else:
            V = [re.sub('[*][0-9][0-9]', '', v) for v in dat_heavy['v_call']]
        J = [re.sub('[*][0-9][0-9]', '', j) for j in dat_heavy['j_call']]
    elif by_alleles:
        if 'v_call_genotyped' in dat_heavy.columns:
            V = [v for v in dat_heavy['v_call_genotyped']]
        else:
            V = [v for v in dat_heavy['v_call']]
        J = [j for j in dat_heavy['j_call']]
    else:
        raise ValueError("by_alleles only accepts boolean values or None, with None defaulting to False.")

    # collapse the alleles to just genes
    V = [','.join(list(set(v.split(',')))) for v in V]
    J = [','.join(list(set(j.split(',')))) for j in J]

    if clustering_by is None or clustering_by is 'aa':
        junction = dict(zip(dat_heavy.index, dat_heavy['junction_aa']))
        if 'junction_aa_length' not in dat_heavy.columns:
            junction_length = [len(str(l)) for l in dat_heavy['junction_aa']]
        else:
            junction_length = [l for l in dat_heavy['junction_aa_length']]
        junction_length_dict = dict(zip(dat_heavy.index, junction_length))
    elif clustering_by == 'nt':
        junction = dict(zip(dat_heavy.index, dat_heavy['junction']))
        if 'junction_length' not in dat_heavy.columns:
            junction_length = [len(str(l)) for l in dat_heavy['junction']]
        else:
            junction_length = [l for l in dat_heavy['junction_length']]
        junction_length_dict = dict(zip(dat_heavy.index, junction_length))
    else:
        raise ValueError("clustering_by only accepts string values 'aa', 'nt' or None, with None defaulting to 'aa'.")
    # Create a dictionary and group sequence ids with same V and J genes
    V_J = dict(zip(dat_heavy.index, zip(V,J)))
    vj_grp = defaultdict(list)
    for key, val in sorted(V_J.items()):
        vj_grp[val].append(key)
    # and now we split the groups based on lengths of the junctions
    vj_len_grp = Tree()
    junction_grp = Tree()
    for g in vj_grp:
        # first, obtain what's the unique lengths
        jlen = []
        for contig_id in vj_grp[g]:
            jlen.append(junction_length_dict[contig_id])
        setjlen = list(set(jlen))
        # then for each unique length, we add the contigs to a new tree if it matches the length
        # and also the actual junction sequences into a another one
        for s in setjlen:
            for contig_id in vj_grp[g]:
                jlen_ = junction_length_dict[contig_id]
                if jlen_ == s:
                    vj_len_grp[g][s][contig_id].value = 1
                    junction_grp[g][s][junction[contig_id]].value = 1
                    for c in [contig_id]:
                        vj_len_grp[g][s][c] = junction[c]
    clones = Tree()
    # for each junction group, calculate the hamming distance matrix
    for g in tqdm(junction_grp, desc = 'Finding clones based on heavy chains '):
        for l in junction_grp[g]:
            junction_ = list(junction_grp[g][l])
            tdarray = np.array(junction_).reshape(-1,1)
            d_mat = squareform(pdist(tdarray, lambda x,y: hamming(x[0],y[0])))
            # then calculate what the acceptable threshold is for each length of sequence
            tr = math.floor(int(l)*(1-identity))
            # convert diagonal and upper triangle to zeroes
            d_mat = np.tril(d_mat)
            np.fill_diagonal(d_mat, 0)
            # get the coordinates/indices of junctions to match against the threshold later
            indices_temp = []
            indices = []
            indices_temp = [list(x) for x in np.tril_indices_from(d_mat)]
            indices = list(zip(indices_temp[0], indices_temp[1]))
            # if there's more than 1 contig, remove the diagonal
            if len(indices) > 1:
                for pairs in indices:
                    # remove diagonals
                    if pairs[0] == pairs[1]:
                        indices.remove(pairs)
            indices_j = []
            # use the coordinates/indices to retrieve the junction sequences
            for p in range(0, len(indices)):
                a1, b1 = indices[p]
                indices_j.append(junction_[a1])
                indices_j.append(junction_[b1])
            # retain only the unique sequences
            indices_j_f = list(set(indices_j))
            # convert the distance matrix to coordinate (source) and distance (target) and create it as a dictionary
            source, target = d_mat.nonzero()
            source_target = list(zip(source.tolist(), target.tolist()))
            if len(source) == 0 & len(target) == 0:
                source_target = list([(0,0)])
            dist = {}
            for st in source_target:
                dist.update({st:d_mat[st]})
            cm1 = []
            cm2 = []
            cm3 = []
            # now to calculate which contigs to group
            tr2 = min(dist.values())
            if tr2 <= tr:
                for d in dist:
                    # if the distance is equal to the minimum distance, which is lesser than or equal to the threshold. I will add the sequence to cm1, meaning they passed and should be a clone.
                    if dist[d] == tr2:
                        cm1.append(d)
                    else:
                        # otherwise add into cm3, which later on will get split into individual clones, despite being part of the acceptable threshold for this group of sequences
                        cm3.append(d)
            else:
                # if the value is greater than the the minimum threshold, I will add the entire group into cm3, so they shouldn't be grouped together as a clone
                for d in dist:
                    cm3.append(d)
            # This is based on a simplest scenario e.g. in 3x3 distance matrix where the acceptable distance threshold is 2.
            # SeqA is 1 sequence different from SeqB
            # SeqB is 2 sequences different from SeqC
            # SeqC can only be 2 sequences different from SeqA otherwise it will violate the pair-wise distance matrix calculation
            # Because the distance between SeqA and SeqB is below the threshold, and is essentially the smallest value, they will be grouped in cm1
            # While SeqC is below the acceptable threshold, really it's more different than SeqA/B 'clone_id'. Hence i'm splitting up here.
            # If there's more than 3 sequences in a particular distance matrix, it gets slightly complicated.
            # so i repeat the calculation to only catch those with more than 3 sequences
            if len(dist) > 3:
                for d in dist:
                    # This time, i will rely on the main threshold value; if the distance lesser than the overall threshold, i will add it to cm2.
                    if dist[d] < tr:
                        cm2.append(d)
                    else:
                        cm3.append(d)
                # any sequences that were already binned above in cm1, it will be removed from here.
                cm2 = list(set(cm2) ^ set(cm1))
                cm3 = list(set(cm3) ^ set(cm1))

            # now to actually retrieve the junction/contigs properly through a series of appending and removing
            j_list1 = []
            if len(cm1) > 0:
                for i in range(0, len(cm1)):
                    a, b = cm1[i]
                    j_list1.append(junction_[a])
                    j_list1.append(junction_[b])
                j_list1 = list(set(j_list1))

            j_list2 = []
            if len(cm3) > 0:
                for i in range(0, len(cm2)):
                    a, b = cm2[i]
                    j_list2.append(junction_[a])
                    j_list2.append(junction_[b])
                j_list2 = list(set(j_list2))

            j_list3 = []
            if len(cm3) > 0:
                for i in range(0, len(cm3)):
                    a, b = cm3[i]
                    j_list3.append(junction_[a])
                    j_list3.append(junction_[b])
                j_list3 = list(set(j_list3))

            # this is to catch the overlapping junctions appearing in the lists because of the sequential appending
            for jl3_1 in j_list1:
                if jl3_1 in j_list3:
                    j_list3.remove(jl3_1)
            for jl2_1 in j_list1:
                if jl2_1 in j_list2:
                    j_list2.remove(jl2_1)
            for jl3_2 in j_list2:
                if jl3_2 in j_list3:
                    j_list3.remove(jl3_2)

            j_list3 = [i.split() for i in list(set(j_list3))]

            if len(j_list1) > 0:
                clones[g][l][str(0)] = j_list1
            if len(j_list2) > 0:
                clones[g][l][str(1)] = j_list2
            if len(j_list3) > 0:
                for c in range(0, len(j_list3)):
                    # the +2 here is so that the numbers come up after 1. It doesn't matter because i will reformat the clone ID numbers later
                    clones[g][l][str(c+2)] = j_list3[c]

    clone_dict = {}
    # now to retrieve the contig ids that are grouped together
    cid = Tree()
    for g in clones:
        for l in clones[g]:
            # retrieve the clone 'numbers'
            for c in clones[g][l]:
                grp_junction = clones[g][l][c]
                for key, value in vj_len_grp[g][l].items():
                    if value in grp_junction:
                        cid[g][l][c][key].value = 1
    # rename clone ids - get dictionaries step by step
    first_key = []
    for k1 in cid.keys():
        first_key.append(k1)
    first_key = list(set(first_key))
    first_key_dict = dict(zip(first_key, range(1,len(first_key)+1)))
    # and now for the middle key
    for g in cid:
        second_key = []
        for k2 in cid[g].keys():
            second_key.append(k2)
        second_key = list(set(second_key))
        second_key_dict = dict(zip(second_key, range(1,len(second_key)+1)))
        for l in cid[g]:
            # and now for the last key
            third_key = []
            for k3 in cid[g][l].keys():
                third_key.append(k3)
            third_key = list(set(third_key))
            third_key_dict = dict(zip(third_key, range(1,len(third_key)+1)))
            for key, value in dict(cid[g][l]).items():
                vL = []
                for v in value:
                    if type(v) is int:
                        break
                    # instead of converting to another tree, i will just make it a dictionary
                    clone_dict[v] = str(first_key_dict[g])+'_'+str(second_key_dict[l])+'_'+str(third_key_dict[key])
    # add it to the original dataframes
    dat_heavy['clone_id'] = pd.Series(clone_dict)
    hclone = dict(zip(dat_heavy['cell_id'], dat_heavy['clone_id']))
    hlclone = dict(zip(dat['sequence_id'], [hclone[c] for c in dat['cell_id']]))
    dat['clone_id'] = pd.Series(hlclone)
    # repeat this process for the light chains within each clone, but only for those with more than 1 light chains in a clone
    dat_light = dat[~(dat['locus'] == 'IGH')]
    # retrieve the J genes and J genes
    for c in tqdm(list(set(dat_light['clone_id'])), desc = 'Refining clone assignment based on light chain pairing '):
        dat_light_c = dat_light[dat_light['clone_id'] == c]
        if dat_light_c.shape[0] > 1:
            if by_alleles is None or ~by_alleles:
                if 'v_call_genotyped' in dat_light_c.columns:
                    Vlight = [re.sub('[*][0-9][0-9]', '', v) for v in dat_light_c['v_call_genotyped']]
                else:
                    Vlight = [re.sub('[*][0-9][0-9]', '', v) for v in dat_light_c['v_call']]
                Jlight = [re.sub('[*][0-9][0-9]', '', j) for j in dat_light_c['j_call']]
            elif by_alleles:
                if 'v_call_genotyped' in dat_light_c.columns:
                    Vlight = [v for v in dat_light_c['v_call_genotyped']]
                else:
                    Vlight = [v for v in dat_light_c['v_call']]
                Jlight = [j for j in dat_light_c['j_call']]
            else:
                raise ValueError("by_alleles only accepts boolean values or None, with None defaulting to False.")
            # collapse the alleles to just genes
            Vlight = [','.join(list(set(v.split(',')))) for v in Vlight]
            Jlight = [','.join(list(set(j.split(',')))) for j in Jlight]
            if clustering_by is None or clustering_by is 'aa':
                junction = dict(zip(dat_light_c.index, dat_light_c['junction_aa']))
                if 'junction_aa_length' not in dat_light.columns:
                    junction_length = [len(str(l)) for l in dat_light_c['junction_aa']]
                else:
                    junction_length = [l for l in dat_light_c['junction_aa_length']]
                junction_length_dict = dict(zip(dat_light_c.index, junction_length))
            elif clustering_by == 'nt':
                junction = dict(zip(dat_light_c.index, dat_light_c['junction']))
                if 'junction_length' not in dat_light_c.columns:
                    junction_length = [len(str(l)) for l in dat_light_c['junction']]
                else:
                    junction_length = [l for l in dat_light_c['junction_length']]
                junction_length_dict = dict(zip(dat_light_c.index, junction_length))
            else:
                raise ValueError("clustering_by only accepts string values 'aa', 'nt' or None, with None defaulting to 'aa'.")

            # Create a dictionary and group sequence ids with same V and J genes
            V_Jlight = dict(zip(dat_light_c.index, zip(Vlight,Jlight)))
            vj_lightgrp = defaultdict(list)
            for key, val in sorted(V_Jlight.items()):
                vj_lightgrp[val].append(key)
            # and now we split the groups based on lengths of the junctions
            vj_len_lightgrp = Tree()
            junction_lightgrp = Tree()
            for g in vj_lightgrp:
                # first, obtain what's the unique lengths
                jlen = []
                for contig_id in vj_lightgrp[g]:
                    jlen.append(junction_length_dict[contig_id])
                setjlen = list(set(jlen))
                # then for each unique length, we add the contigs to a new tree if it matches the length
                # and also the actual junction sequences into a another one
                for s in setjlen:
                    for contig_id in vj_lightgrp[g]:
                        jlen_ = junction_length_dict[contig_id]
                        if jlen_ == s:
                            vj_len_lightgrp[g][s][contig_id].value = 1
                            junction_lightgrp[g][s][junction[contig_id]].value = 1
                            for c in [contig_id]:
                                vj_len_lightgrp[g][s][c] = junction[c]

            clones_light = Tree()
            for g in junction_lightgrp:
                for l in junction_lightgrp[g]:
                    junction_ = list(junction_lightgrp[g][l])
                    tdarray = np.array(junction_).reshape(-1,1)
                    d_mat = squareform(pdist(tdarray, lambda x,y: hamming(x[0],y[0])))
                    tr = math.floor(int(l)*(1-identity))
                    d_mat = np.tril(d_mat)
                    np.fill_diagonal(d_mat, 0)
                    indices_temp = []
                    indices = []
                    indices_temp = [list(x) for x in np.tril_indices_from(d_mat)]
                    indices = list(zip(indices_temp[0], indices_temp[1]))
                    if len(indices) > 1:
                        for pairs in indices:
                            if pairs[0] == pairs[1]:
                                indices.remove(pairs)
                    indices_j = []
                    for p in range(0, len(indices)):
                        a1, b1 = indices[p]
                        indices_j.append(junction_[a1])
                        indices_j.append(junction_[b1])
                    indices_j_f = list(set(indices_j))
                    source, target = d_mat.nonzero()
                    source_target = list(zip(source.tolist(), target.tolist()))
                    if len(source) == 0 & len(target) == 0:
                        source_target = list([(0,0)])
                    dist = {}
                    for st in source_target:
                        dist.update({st:d_mat[st]})
                    cm1 = []
                    cm2 = []
                    cm3 = []
                    tr2 = min(dist.values())
                    if tr2 <= tr:
                        for d in dist:
                            if dist[d] == tr2:
                                cm1.append(d)
                            else:
                                cm3.append(d)
                    else:
                        for d in dist:
                            cm3.append(d)

                    if len(dist) > 3:
                        for d in dist:
                            if dist[d] < tr:
                                cm2.append(d)
                            else:
                                cm3.append(d)
                        cm2 = list(set(cm2) ^ set(cm1))
                        cm3 = list(set(cm3) ^ set(cm1))
                    j_list1 = []
                    if len(cm1) > 0:
                        for i in range(0, len(cm1)):
                            a, b = cm1[i]
                            j_list1.append(junction_[a])
                            j_list1.append(junction_[b])
                        j_list1 = list(set(j_list1))
                    j_list2 = []
                    if len(cm3) > 0:
                        for i in range(0, len(cm2)):
                            a, b = cm2[i]
                            j_list2.append(junction_[a])
                            j_list2.append(junction_[b])
                        j_list2 = list(set(j_list2))
                    j_list3 = []
                    if len(cm3) > 0:
                        for i in range(0, len(cm3)):
                            a, b = cm3[i]
                            j_list3.append(junction_[a])
                            j_list3.append(junction_[b])
                        j_list3 = list(set(j_list3))
                    for jl3_1 in j_list1:
                        if jl3_1 in j_list3:
                            j_list3.remove(jl3_1)
                    for jl2_1 in j_list1:
                        if jl2_1 in j_list2:
                            j_list2.remove(jl2_1)
                    for jl3_2 in j_list2:
                        if jl3_2 in j_list3:
                            j_list3.remove(jl3_2)
                    j_list3 = [i.split() for i in list(set(j_list3))]
                    if len(j_list1) > 0:
                        clones_light[g][l][str(0)] = j_list1
                    if len(j_list2) > 0:
                        clones_light[g][l][str(1)] = j_list2
                    if len(j_list3) > 0:
                        for c in range(0, len(j_list3)):
                            clones_light[g][l][str(c+2)] = j_list3[c]
            clone_dict_light = {}
            cid_light = Tree()
            for g in clones_light:
                for l in clones_light[g]:
                    for c in clones_light[g][l]:
                        grp_junction = clones_light[g][l][c]
                        for key, value in vj_len_lightgrp[g][l].items():
                            if value in grp_junction:
                                cid_light[g][l][c][key].value = 1
            first_key = []
            for k1 in cid_light.keys():
                first_key.append(k1)
            first_key = list(set(first_key))
            first_key_dict = dict(zip(first_key, range(1,len(first_key)+1)))
            for g in cid_light:
                second_key = []
                for k2 in cid_light[g].keys():
                    second_key.append(k2)
                second_key = list(set(second_key))
                second_key_dict = dict(zip(second_key, range(1,len(second_key)+1)))
                for l in cid_light[g]:
                    third_key = []
                    for k3 in cid_light[g][l].keys():
                        third_key.append(k3)
                    third_key = list(set(third_key))
                    third_key_dict = dict(zip(third_key, range(1,len(third_key)+1)))
                    for key, value in dict(cid_light[g][l]).items():
                        vL = []
                        for v in value:
                            if type(v) is int:
                                break
                            clone_dict_light[v] = str(first_key_dict[g])+'_'+str(second_key_dict[l])+'_'+str(third_key_dict[key])
            lclones = list(clone_dict_light.values())
            # will just update the main dat directly
            if len(list(set(lclones))) > 1:
                lclones_dict = dict(zip(sorted(list(set(lclones))), [str(x) for x in range(1,len(list(set(lclones)))+1)]))
                renamed_clone_dict_light = {}
                for key, value in clone_dict_light.items():
                    renamed_clone_dict_light[key] = lclones_dict[value]
                dat.loc[renamed_clone_dict_light.keys(), 'clone_id'] = dat.loc[renamed_clone_dict_light.keys(), 'clone_id'] + '_' + pd.Series(renamed_clone_dict_light)

    if os.path.isfile(str(self)):
        dat.to_csv("{}/{}_clone.tsv".format(os.path.dirname(self), os.path.basename(self).split('.tsv')[0]), sep = '\t', index = False)
    else:
        if write_out:
            if outdir is None:
                outDir = 'dandelion/data'
            else:
                if outdir.endswith('/'):
                    outDir = str(outdir).strip('/')
                else:
                    outDir = str(outdir)

            if not os.path.exist(outDir):
                os.makedirs(outDir)
            if outFilePrefix is not None:
                dat.to_csv("{}/{}_clone.tsv".format(outDir, str(outFilePrefix)), sep = '\t', index = None)
            elif outFilePrefix is None:
                dat.to_csv("{}/{}_clone.tsv".format(outDir, 'dandelion_find_clones'), sep = '\t', index = None)
        else:
            pass

    if self.__class__ == Dandelion:
        self.__init__(data = dat)
    else:
        out = Dandelion(data = dat)

    logg.info(' finished', time=start,
        deep=('added to Dandelion class object: \n'
        '   \'data\', contig-indexed clone table\n'))

        return(out)

def generate_network(self, distance_mode='weighted', aa_or_nt=None, clones_sep = None, weights = None, layout_option = None, *args, **kwds):
    """
    Generates a levenshtein distance network based on gapped full length sequences for heavy and light chain(s). 
    The distance matrices are then combined into a singular matrix where a minimum spanning tree will be constructed per clone group specified by separator in `clones_sep` option.

    Parameters
    ----------
    data : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones have been determined.
    distance_mode : str
        The mode of calculating joint distance matrix for heavy and light chains. Default is 'weighted'. If 'simple', a simple sum operation will be used. If 'weighted', depending on whether `weights` option is provided, it will scale each layer to range of 0..1 to bring the multiple layers of data into a single analysis.
    aa_or_nt : str, optional
        Option accepts 'aa', 'nt' or None, with None defaulting to 'aa'. Determines whether amino acid or nucleotide sequences will be used for calculating distances.
    clones_sep: tuple(int, str)
        A tuple containing how the clone groups should be extracted. None defaults to (0, '_').
    weights : tuple, optional
        A tuple containing weights to scale each layer. default is None where each layer is scaled evenly i.e. 1/number of layers.
    layout_option : str, optional
        choice of layout algorithm. None defaults to fruchterman reingold layout.
    *args and **kwds
        passed to `igraph.graph.layout <https://igraph.org/python/doc/igraph.Graph-class.html>`__.
    Returns
    ----------
        `Dandelion` object with `.distance`, `.edges`, `.layout`, `.graph` initialized.
    """
    start = logg.info('Generating network')
    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    else:
        dat = load_data(self)
    if 'clone_id' not in dat.columns:
        raise TypeError('Data does not contain clone information. Please run find_clones.')

    # re-initiate a Dandelion class object
    out = Dandelion(dat)

    # calculate distance
    dat_h = dat[dat['locus'] == 'IGH']
    dat_l = dat[dat['locus'].isin(['IGK', 'IGL'])]
    if aa_or_nt is None or aa_or_nt is 'aa':
        seq_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h['sequence_alignment_aa'])))
        seq_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l['sequence_alignment_aa'])))
    elif aa_or_nt == 'nt':
        seq_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h['sequence_alignment'])))
        seq_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l['sequence_alignment'])))
    else:
        raise ValueError("aa_or_nt only accepts string values 'aa', 'nt' or None, with None defaulting to 'aa'.")

    # So first, create a data frame to hold all possible (full) sequences split by heavy (only 1 possible) and light (multiple possible)
    dat_seq = pd.DataFrame.from_dict(seq_h, orient = 'index', columns = ['cell_id', 'heavy'])
    dat_seq.set_index('cell_id', inplace = True)
    light_seq_tree = Tree()
    for key, value in seq_l.items():
        k, v = value
        light_seq_tree[k][key] = v
    light_seq_tree2 = Tree()
    for g in light_seq_tree:
        second_key = []
        for k2 in light_seq_tree[g].keys():
            second_key.append(k2)
        second_key = list(set(second_key))
        second_key_dict = dict(zip(second_key, range(0,len(second_key))))
        for key, value in light_seq_tree[g].items():
            light_seq_tree2[g][second_key_dict[key]] = value
    dat_seq['light'] = pd.Series(light_seq_tree2)
    tmp_dat = dat_seq['light'].apply(pd.Series)
    tmp_dat.columns = ['light_' + str(c) for c in tmp_dat.columns]
    dat_seq = dat_seq.merge(tmp_dat, left_index = True, right_index = True)
    dat_seq = dat_seq[['heavy'] + [str(c) for c in tmp_dat.columns]]

    # calculate a distance matrix for all vs all and this can be referenced later on to extract the distance between the right pairs
    dmat = Tree()
    for x in tqdm(dat_seq.columns, desc = 'Calculating distances... '):
        seq_list = []
        seq_list = [y for y in dat_seq[x]]
        tdarray = np.array(seq_list).reshape(-1,1)
        d_mat = squareform(pdist(tdarray,lambda x,y: Levenshtein.distance(x[0],y[0])))
        dmat[x] = d_mat
    dist_mat_list = [dmat[x] for x in dmat if type(dmat[x]) is np.ndarray]
    
    n_ = len(dist_mat_list)
    if distance_mode == 'simple':
        total_dist = np.sum(dist_mat_list,axis=0)
    if distance_mode == 'weighted':
        weighted_matrix = []
        if weights is None:
            for w in range(0, n_):
                weighted_matrix.append(1/n_ * dist_mat_list[w])
            total_dist = sum(weighted_matrix)
        else:
            if len(weights) == n_:                
                for w in range(0, n_):
                    weighted_matrix.append(weights[w] * dist_mat_list[w])
                total_dist = sum(weighted_matrix)
            else:
                raise IndexError('Length of provided weights should be %s.' % int(n_))

    # generate edge list
    tmp_totaldist = pd.DataFrame(total_dist, index = out.metadata.index, columns = out.metadata.index)
    tmp_clusterdist = Tree()
    for i in out.metadata.index:
        cx = out.metadata.loc[i,'clone_group_id']
        tmp_clusterdist[cx][i].value = 1
    tmp_clusterdist2 = {}
    for x in tmp_clusterdist:
        tmp_clusterdist2[x] = list(tmp_clusterdist[x])
    cluster_dist = {}
    for x in tmp_clusterdist2:
        dist_mat_ = tmp_totaldist.loc[tmp_clusterdist2[x], tmp_clusterdist2[x]]
        s1, s2 = dist_mat_.shape
        if s1 > 1 and s2 >1:
            cluster_dist[x] = dist_mat_
    # to improve the visulisation and plotting efficiency, i will build a minimum spanning tree for each group/clone to connect the shortest path
    mst_tree = mst(cluster_dist)
    sleep(0.5)

    edge_list = Tree()
    for c in tqdm(mst_tree, desc = 'Generating edge list '):
        G = nx.from_pandas_adjacency(mst_tree[c], create_using=nx.MultiDiGraph())
        G.edges(data=True)
        edge_list[c] = nx.to_pandas_edgelist(G)
    sleep(0.5)
    clone_ref = dict(out.metadata['clone_id'])
    tmp_clone_tree = Tree()
    for x in out.metadata.index:
        tmp_clone_tree[clone_ref[x]][x].value = 1
    tmp_clone_tree2 = Tree()
    for x in tmp_clone_tree:
        tmp_clone_tree2[x] = list(tmp_clone_tree[x])

    tmp_clone_tree3 = Tree()
    for x in tmp_clone_tree2:
        tmp_ = pd.DataFrame(index = tmp_clone_tree2[x], columns = tmp_clone_tree2[x])
        tmp_ = pd.DataFrame(np.tril(tmp_) + 1, index = tmp_clone_tree2[x], columns = tmp_clone_tree2[x])
        tmp_.fillna(0, inplace = True)
        tmp_clone_tree3[x] = tmp_

    # here I'm using a temporary edge list to catch all cells that were identified as clones to forecfully link them up if they were clipped off during the mst step
    tmp_edge_list = Tree()
    for c in tqdm(tmp_clone_tree3, desc = 'Linking edges '):
        G = nx.from_pandas_adjacency(tmp_clone_tree3[c], create_using=nx.MultiDiGraph())
        G.edges(data=True)
        tmp_edge_list[c] = nx.to_pandas_edgelist(G)

    edge_listx = pd.concat([edge_list[x] for x in edge_list])
    edge_listx.index = [(s, t) for s, t in zip(edge_listx['source'],edge_listx['target'])]

    tmp_edge_listx = pd.concat([tmp_edge_list[x] for x in tmp_edge_list])
    tmp_edge_listx.index = [(s, t) for s, t in zip(tmp_edge_listx['source'], tmp_edge_listx['target'])]

    edge_list_final = edge_listx.combine_first(tmp_edge_listx)

    for idx in edge_list_final.index:
        edge_list_final.loc[idx, 'weight'] = tmp_totaldist.loc[idx[0], idx[1]]
    # return the edge list
    edge_list_final.reset_index(drop = True, inplace = True)

    # and finally the vertex list which is super easy
    vertices = pd.DataFrame(out.metadata.index)

    # and now to actually generate the network
    edges = [tuple(e) for e in edge_list_final.values]
    edges_network = [tuple((e[0], e[1])) for e in edges]
    graph = igraph.Graph.Formula()
    graph.add_vertices(vertices['cell_id'])
    graph.add_edges(edges_network) # may take a very long time if there's too many

    if layout_option is not None:
        layout = graph.layout_fruchterman_reingold()
    else:
        layout = graph.layout(layout_option, *args, **kwds)
    for x in out.metadata.columns:
        graph.vs[x] = out.metadata[x]
    graph.es['width'] = [0.8/(int(e[2]) + 1) for e in edges]

    if self.__class__ == Dandelion:
        self.__init__(data = dat, distance = dmat, edges = edge_list_final, layout = layout, graph = graph)
    else:
        out = Dandelion(data = dat, distance = dmat, edges = edge_list_final, layout = layout, graph = graph)
    
    logg.info(' finished', time=start,
        deep=('added to Dandelion class object: \n'
        '   \'data\', contig-indexed clone table\n'
        '   \'metadata\', cell-indexed clone table\n'
        '   \'distance\', heavy and light chain distance matrices\n'
        '   \'edges\', network edges\n'
        '   \'layout\', network layout\n'
        '   \'graph\', network'))
    
    return(out)

def mst(mat):
    """
    Construct minimum spanning tree based on supplied matrix in dictionary.

    Parameters
    ----------
    mat : dict
        Dictionary containing numpy ndarrays.
    Returns
    ----------        
        Dandelion `Tree` object holding DataFrames of constructed minimum spanning trees.
    """
    mst_tree = Tree()
    for c in mat:
        mst_tree[c] = pd.DataFrame(minimum_spanning_tree(np.triu(mat[c])).toarray().astype(int), index = mat[c].index, columns = mat[c].columns)
    return(mst_tree)

def transfer_network(self, network, keep_raw = True, neighbors_key = None):
    """
    Transfer data in `Dandelion` slots to `AnnData` object, updating the `.obs`, `.uns`, `.obsm` and `.raw` slots with metadata and network.

    Parameters
    ----------
    self : AnnData
        `AnnData` object
    network : Dandelion
        `Dandelion` object
    keep_raw : bool
        If True, will transfer the existing `.uns` slot to `.raw.uns`.
    neighbors_key : str, optional
        key for 'neighbors' slot in `.uns`
    Returns
    ----------
        `AnnData` object with updated `.obs` `.uns`, `.obsm` (and `.raw`) slots with data from `Dandelion` object.

    """
    start = logg.info('Transferring network')
    G = nx.from_pandas_edgelist(network.edges, create_using=nx.MultiDiGraph(), edge_attr='weight')
    distances = nx.to_pandas_adjacency(G, dtype = np.float32, weight='weight')
    connectivities = nx.to_pandas_adjacency(G, dtype = np.float32, weight=None)
    df_connectivities = pd.DataFrame(index = self.obs.index, columns = self.obs.index)
    df_distances = pd.DataFrame(index = self.obs.index, columns = self.obs.index)
    for x in connectivities.columns:
        df_connectivities[x] = pd.Series(connectivities[x])
    for x in distances.columns:
        df_distances[x] = pd.Series(distances[x])
    for x in df_distances.columns:
        df_distances[x] = df_distances[x].apply(lambda x: 5/(x + 1))
    df_connectivities.fillna(0, inplace = True)
    df_distances.fillna(0, inplace = True)
    df_connectivities_ = scipy.sparse.csr_matrix(df_connectivities.values, dtype = np.float32)
    df_distances_ = scipy.sparse.csr_matrix(df_distances.values, dtype = np.float32)

    if neighbors_key is None:
        neighbors_key = "neighbors"
    if neighbors_key not in self.uns:
        raise ValueError("`edges=True` requires `pp.neighbors` to be run before.")
    if keep_raw:
        self.raw.uns = copy.deepcopy(self.uns)
        self.uns[neighbors_key]['connectivities'] = df_connectivities_
        self.uns[neighbors_key]['distances'] = df_distances_
        self.uns[neighbors_key]['params'] = {'method':'bcr'}
    else:
        self.uns[neighbors_key]['connectivities'] = df_connectivities_
        self.uns[neighbors_key]['distances'] = df_distances_
        self.uns[neighbors_key]['params'] = {'method':'bcr'}

    for x in network.metadata.columns:
        self.obs[x] = pd.Series(network.metadata[x])

    tmp = self.obs.copy()
    coord = pd.DataFrame(np.array(network.layout), index = network.metadata.index)
    for x in coord.columns:
        tmp[x] = coord[x]
    tmp[[1]] = tmp[[1]]*-1
    X_bcr = np.array(tmp[[0,1]], dtype = np.float32)
    self.obsm['X_bcr'] = X_bcr
    if keep_raw:
        logg.info(' finished', time=start,
            deep=('added to `.uns[\''+neighbors_key+'\']`\n'
            '   \'distances\', cluster-weighted adjacency matrix\n'
            '   \'connectivities\', cluster-weighted adjacency matrix\n'
            'stored original .uns in .raw'))
    else:
        logg.info(' finished', time=start,
            deep=('added to `.uns[\''+neighbors_key+'\']`\n'
            '   \'distances\', cluster-weighted adjacency matrix\n'
            '   \'connectivities\', cluster-weighted adjacency matrix'))

def quantify_mutations(self, split_locus = False, region_definition=None, mutation_definition=None, frequency=True, combine=True):
    """
    Runs basic mutation load analysis implemented in `shazam <https://shazam.readthedocs.io/en/stable/vignettes/Mutation-Vignette/>`__. 

    Parameters
    ----------
    self : Dandelion
        `Dandelion` object
    split_locus : bool
        whether to return the results for heavy chain and light chain separately. Default is False.
    region_definition : str, optional
        passed to shazam's `observedMutations <https://shazam.readthedocs.io/en/stable/topics/observedMutations/>`__.
    mutation_definition : str, optional
        passed to shazam's `observedMutations <https://shazam.readthedocs.io/en/stable/topics/observedMutations/>`__.
    frequency
        whether to return the results a frequency or counts. Default is True (frequency).
    combine
        whether to return the results for replacement and silent mutations separately (False). Default is True (sum).
    Returns
    ----------
        `Dandelion` object with updated `.metadata` slot.
    """
    sh = importr('shazam')
    dat = load_data(self.data)
    warnings.filterwarnings("ignore")

    if region_definition is None:
        reg_d = NULL
    else:
        reg_d = data(sh).fetch(region_definition)

    if mutation_definition is None:
        mut_d = NULL
    else:
        mut_d = data(sh).fetch(mutation_definition)

    if split_locus is False:
        try:
            dat_r = pandas2ri.py2rpy(dat)
        except:
            dat = dat.fillna('')
            dat_r = pandas2ri.py2rpy(dat)

        results = sh.observedMutations(dat_r, sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask", regionDefinition = reg_d, mutationDefinition = mut_d, frequency = frequency, combine = combine)
        pd_df = pandas2ri.rpy2py_dataframe(results)
    else:
        dat_h = dat[dat['locus'] == 'IGH']
        dat_l = dat[dat['locus'].isin(['IGK', 'IGL'])]

        try:
            dat_h_r = pandas2ri.py2rpy(dat_h)
        except:
            dat_h = dat_h.fillna('')
            dat_h_r = pandas2ri.py2rpy(dat_h)

        try:
            dat_l_r = pandas2ri.py2rpy(dat_l)
        except:
            dat_l = dat_l.fillna('')
            dat_l_r = pandas2ri.py2rpy(dat_l)

        results_h = sh.observedMutations(dat_h_r, sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask", regionDefinition = reg_d, mutationDefinition = mut_d, frequency = frequency, combine = combine)
        results_l = sh.observedMutations(dat_l_r, sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask", regionDefinition = reg_d, mutationDefinition = mut_d, frequency = frequency, combine = combine)
        pd_df_h = pandas2ri.rpy2py_dataframe(results_h)
        pd_df_l = pandas2ri.rpy2py_dataframe(results_l)
        pd_df = pd.concat([pd_df_h, pd_df_l])

    pd_df.set_index('sequence_id', inplace = True, drop = False)
    cols_to_return = pd_df.columns.difference(dat.columns) # this doesn't actually catch overwritten columns
    if len(cols_to_return) < 1:
        cols_to_return = list(filter(re.compile("mu_.*").match, [c for c in pd_df.columns]))
    else:
        cols_to_return = cols_to_return
    res = {}
    for x in cols_to_return:
        res[x] = list(pd_df[x])
        self.data[x] = [str(r) for r in res[x]] # TODO: str will make it work for the back and forth conversion with rpy2. but maybe can use a better option?
    if split_locus is False:
        metadata_ = self.data[['cell_id']+list(cols_to_return)]
    else:
        metadata_ = self.data[['locus', 'cell_id']+list(cols_to_return)]
    
    for x in cols_to_return:
        metadata_[x] = metadata_[x].astype(np.float32)

    if split_locus is False:
        metadata_ = metadata_.groupby('cell_id').sum()
    else:
        metadata_ = metadata_.groupby(['locus','cell_id']).sum()
        metadatas = []
        for x in list(set(self.data['locus'])):
            tmp = metadata_.iloc[metadata_.index.isin([x], level='locus'),:]
            tmp.index = tmp.index.droplevel()
            tmp.columns = [c+'_'+str(x) for c in tmp.columns]
            metadatas.append(tmp)
        metadata_ = functools.reduce(lambda x, y: pd.merge(x, y, left_index = True, right_index = True, how = 'outer'), metadatas)

    metadata_.index.name = None
    if self.metadata is None:
        self.metadata = metadata_
    else:
        for x in metadata_.columns:
            self.metadata[x] = pd.Series(metadata_[x])

def calculate_threshold(self, manual_threshold=None, model=None, normalize_method=None, threshold_method=None, edge=None, cross=None, subsample=None, threshold_model=None, cutoff=None, sensitivity=None, specificity=None, ncpu=None, plot=True, plot_group=None,  figsize=(4.5, 2.5), *args):
    """
    Calculating nearest neighbor distances for tuning clonal assignment with `shazam <https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/>`__.
    
    Runs the following:        
    distToNearest
        Get non-zero distance of every heavy chain (IGH) sequence (as defined by sequenceColumn) to its nearest sequence in a partition of heavy chains sharing the same V gene, J gene, and junction length (VJL), or in a partition of single cells with heavy chains sharing the same heavy chain VJL combination, or of single cells with heavy and light chains sharing the same heavy chain VJL and light chain VJL combinations.
    findThreshold
        automtically determines an optimal threshold for clonal assignment of Ig sequences using a vector of nearest neighbor distances. It provides two alternative methods using either a Gamma/Gaussian Mixture Model fit (threshold_method="gmm") or kernel density fit (threshold_method="density").
    
    Parameters
    ----------
    self : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones have been determined.
    manual_threshold : float, optional
        value to manually plot in histogram.
    model : str, optional
        underlying SHM model, which must be one of c("ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "hs1f_compat", "m1n_compat").
    normalize_method : str, optional
        method of normalization. The default is "len", which divides the distance by the length of the sequence group. If "none" then no normalization if performed.
    threshold_method : str, optional
        string defining the method to use for determining the optimal threshold. One of "gmm" or "density".
    edge : float, optional
        upper range as a fraction of the data density to rule initialization of Gaussian fit parameters. Default value is 0.9 (or 90). Applies only when threshold_method="density".
    cross : list, array, optional
        supplementary nearest neighbor distance vector output from distToNearest for initialization of the Gaussian fit parameters. Applies only when method="gmm".
    subsample : int, optional
        maximum number of distances to subsample to before threshold detection.
    threshold_model : str, optional
        allows the user to choose among four possible combinations of fitting curves: "norm-norm", "norm-gamma", "gamma-norm", and "gamma-gamma". Applies only when method="gmm".
    cutoff : str, optional
        method to use for threshold selection: the optimal threshold "opt", the intersection point of the two fitted curves "intersect", or a value defined by user for one of the sensitivity or specificity "user". Applies only when method="gmm".
    sensitivity : float, optional
        sensitivity required. Applies only when method="gmm" and cutoff="user".
    specificity : float, optional
        specificity required. Applies only when method="gmm" and cutoff="user".
    ncpu : int, optional
        number of cpus for parallelization. Default is all available cpus.
    plot : bool
        whether or not to return plot.
    plot_group : str, optional
        determines the fill color and facets.
    figsize : tuple
        size of plot. Default is (4.5, 2.5).
    *args
        passed to shazam's `distToNearest <https://shazam.readthedocs.io/en/stable/topics/distToNearest/>`__.
    Returns
    ----------
        plotnine plot showing histogram of length normalized ham model distance threshold.
    """

    sh = importr('shazam')
    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    elif self.__class__ == pd.DataFrame or os.path.isfile(str(self)):
        dat = load_data(self)
    warnings.filterwarnings("ignore")

    if 'v_call_genotyped' in dat.columns:
        v_call = 'v_call_genotyped'
    else:
        v_call = 'v_call'
    
    if model is None:
        model_ = 'ham'
    else:
        model_ = model
    
    if normalize_method is None:
        norm_ = 'len'
    else:
        norm_ = normalize_method
    
    if threshold_method is None:
        threshold_method_ = "density"
    else:
        threshold_method_ = threshold_method
    
    if subsample is None:
        subsample_ = NULL
    else:
        subsample_ = subsample
    
    if ncpu is None:
        ncpu_ = multiprocessing.cpu_count()
    else:
        ncpu_ = ncpu

    dat_h = dat[dat['locus'] == 'IGH']
    
    try:
        dat_h_r = pandas2ri.py2rpy(dat_h)
    except:
        dat_h = dat_h.fillna('')
        dat_h_r = pandas2ri.py2rpy(dat_h)

    dist_ham = sh.distToNearest(dat_h_r, vCallColumn=v_call, model=model_, normalize=norm_, nproc=ncpu_, *args)
    # Find threshold using density method
    c = rpy2.robjects.StrVector(['dist_nearest'])

    if threshold_method_ is 'density':
        if edge is None:
            edge_ = 0.9
        else:
            edge_ = edge
            
        dist_threshold = sh.findThreshold(dist_ham.rx(True, c), method=threshold_method_, subsample = subsample_, edge = edge_)
    else:
        if threshold_model is None:
            threshold_model_ = "gamma-gamma"
        else:
            threshold_model_ = threshold_model
        
        if cross is None:
            cross_ = NULL
        else:
            cross_ = cross
        
        if cutoff is None:
            cutoff_ = 'optimal'
        else:
            cutoff_ = cutoff
        
        if sensitivity is None:
            sen_ = NULL
        else:
            sen_ = sensitivity
        
        if specificity is None:
            spc_ = NULL
        else:
            spc_ = specificity        
        dist_threshold = sh.findThreshold(dist_ham.rx(True, c), method=threshold_method, model = threshold_model_, cross = cross_, subsample = subsample_, cutoff = cutoff_, sen = sen_, spc = spc_)        

    threshold=np.array(dist_threshold.slots['threshold'])[0]
    
    dist_ham = pandas2ri.rpy2py_dataframe(dist_ham)
    
    if plot:
        options.figure_size = figsize
        if plot_group is None:
            plot_group = 'sample_id'
        else:
            plot_group = plot_group
        if manual_threshold is None:
            tr = threshold
        else:
            tr = manual_threshold            
        p = (ggplot(dist_ham, aes('dist_nearest', fill=str(plot_group)))
             + theme_bw() 
             + xlab("Grouped Hamming distance")
             + ylab("Count")
             + geom_histogram(binwidth = 0.01)
             + geom_vline(xintercept = tr, linetype = "dashed", color="blue", size=0.5)
             + annotate('text', x=tr+0.02, y = 10, label='Threshold:\n' + str(np.around(tr, decimals=2)), size = 8, color = 'Blue', hjust = 'left')
             + facet_grid('~'+str(plot_group), scales="free_y")
             + theme(legend_position = 'none'))        
        return(p)            
    else:
        print('Automatic Threshold : '+str(np.around(threshold, decimals=2), '\n method = '+str(threshold_method)))

def define_clones(self, dist, action = 'set', model = 'ham', norm = 'len', doublets='drop', fileformat='airr', ncpu = None, dirs = None, outFilePrefix = None, verbose = False):
    """
    Find clones using changeo's `DefineClones.py <https://changeo.readthedocs.io/en/stable/tools/DefineClones.html>`__.
    
    Parameters
    ----------
    self : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones have been determined.
    dist : float
        The distance threshold for clonal grouping.
    action : str
        Specifies how to handle multiple V(D)J assignments for initial grouping. Default is 'set'. The “first” action will use only the first gene listed. The “set” action will use all gene assignments and construct a larger gene grouping composed of any sequences sharing an assignment or linked to another sequence by a common assignment (similar to single-linkage).
    model : str
        Specifies which substitution model to use for calculating distance between sequences. Default is 'ham'. The “ham” model is nucleotide Hamming distance and “aa” is amino acid Hamming distance. The “hh_s1f” and “hh_s5f” models are human specific single nucleotide and 5-mer content models, respectively, from Yaari et al, 2013. The “mk_rs1nf” and “mk_rs5nf” models are mouse specific single nucleotide and 5-mer content models, respectively, from Cui et al, 2016. The “m1n_compat” and “hs1f_compat” models are deprecated models provided backwards compatibility with the “m1n” and “hs1f” models in Change-O v0.3.3 and SHazaM v0.1.4. Both 5-mer models should be considered experimental.
    norm : str
        Specifies how to normalize distances. Default is 'len'. 'none' (do not normalize), 'len' (normalize by length), or 'mut' (normalize by number of mutations between sequences).
    doublets : str
        Option to control behaviour when dealing with heavy chain 'doublets'. Default is 'drop'. 'drop' will filter out the doublets while 'count' will retain only the highest umi count contig.
    fileformat : str
        format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.
    ncpu : int, optional
        number of cpus for parallelization. Default is all available cpus.
    dirs : str, optional
        If specified, out file will be in this location.
    outFilePrefix : str, optional
        If specified, the out file name will have this prefix. None defaults to 'dandelion_define_clones'
    verbose : bool
        Whether or not to print the command used in terminal to call DefineClones.py. Default is False.
    Returns
    ----------
        `Dandelion` object with clone_id annotated in `.data` slot and `.metadata` initialized.
    """
    if ncpu is None:
        nproc=multiprocessing.cpu_count()
    else:
        nproc=ncpu

    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    else:
        dat = load_data(self)
    if os.path.isfile(str(self)):
        dat = load_data(self)
    dat_h = dat[dat['locus'] == 'IGH']
    dat_l = dat[dat['locus'].isin(['IGK', 'IGL'])]
    
    if os.path.isfile(str(self)):
        if dirs is None:
            tmpFolder = "{}/tmp".format(os.path.dirname(self))
            outFolder = "{}".format(os.path.dirname(self))
        else:
            tmpFolder = str(dirs).strip('/')+'/tmp'
            outFolder = str(dirs).strip('/')
    else:
        if dirs is None:
            tmpFolder = "dandelion/data/tmp"
            outFolder = "dandelion/data"
        else:
            tmpFolder = str(dirs).strip('/')+'/tmp'
            outFolder = str(dirs).strip('/')

    if not os.path.exists(tmpFolder):
        os.makedirs(tmpFolder)
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)

    if os.path.isfile(str(self)):
        h_file1 = "{}/{}_heavy-clone.tsv".format(tmpFolder, os.path.basename(self).split('.tsv')[0])
        h_file2 = "{}/{}_heavy-clone.tsv".format(outFolder, os.path.basename(self).split('.tsv')[0])
        l_file = "{}/{}_light.tsv".format(tmpFolder, os.path.basename(self).split('.tsv')[0])
        outfile = "{}/{}_clone.tsv".format(outFolder, os.path.basename(self).split('.tsv')[0])
    else:        
        if outFilePrefix is not None:
            out_FilePrefix = outFilePrefix
        else:
            out_FilePrefix = 'dandelion_define_clones'
        h_file1 = "{}/{}_heavy-clone.tsv".format(tmpFolder, out_FilePrefix)
        h_file2 = "{}/{}_heavy-clone.tsv".format(outFolder, out_FilePrefix)
        l_file = "{}/{}_light.tsv".format(tmpFolder, out_FilePrefix)
        outfile = "{}/{}_clone.tsv".format(outFolder, out_FilePrefix)

    dat_h.to_csv(h_file1, sep = '\t', index = False)
    dat_l.to_csv(l_file, sep = '\t', index = False)
    
    if 'germline_alignment_d_mask' not in dat.columns:
        raise ValueError("Missing 'germline_alignment_d_mask' column in input file. Run create_germlines first.")

    if 'v_call_genotyped' in dat.columns:
        v_field = 'v_call_genotyped'
    else:
        v_field = 'v_call'

    cmd = ['DefineClones.py',
            '-d', h_file1,
            '-o', h_file2,
            '--act', action,
            '--model', model,
            '--norm', norm,
            '--dist', str(dist),
            '--nproc', str(nproc),
            '--vf', v_field]
    
    def clusterLinkage(cell_series, group_series):
        """
        Returns a dictionary of {cell_id : cluster_id} that identifies clusters of cells by analyzing their shared
        features (group_series) using single linkage. 
    
        Arguments:
        cell_series (iter): iter of cell ids.
        group_series (iter): iter of group ids.
    
        Returns:
        dict:  dictionary of {cell_id : cluster_id}.
        """
    
        # assign initial clusters
        # initial_dict = {cluster1: [cell1], cluster2: [cell1]}
        initial_dict = {}
        for cell, group in zip(cell_series, group_series):
            try:    
                initial_dict[group].append(cell)
            except KeyError:
                initial_dict[group] = [cell]
                   
        # naive single linkage clustering (ON^2 best case, ON^3 worst case) ...ie for cells with multiple light chains
        # cluster_dict = {cluster1: [cell1, cell2]}, 2 cells belong in same group if they share 1 light chain 
        while True:
            cluster_dict = {}
            for i, group in enumerate(initial_dict.keys()):
                cluster_dict[i] = initial_dict[group]
                for cluster in cluster_dict:
                    # if initial_dict[group] and cluster_dict[cluster] share common cells, add initial_dict[group] to cluster
                    if cluster != i and any(cell in initial_dict[group] for cell in cluster_dict[cluster]):
                        cluster_dict[cluster] = cluster_dict[cluster] + initial_dict[group]
                        del cluster_dict[i]
                        break
            # break if clusters stop changing, otherwise restart 
            if len(cluster_dict.keys()) == len(initial_dict.keys()):
                break
            else:
                initial_dict = cluster_dict.copy()
        
        # invert cluster_dict for return
        assign_dict = {cell:k for k,v in cluster_dict.items() for cell in set(v)}
        
        return assign_dict

    def _lightCluster(heavy_file, light_file, out_file, doublets, fileformat):
        """
        Split heavy chain clones based on light chains
    
        Arguments:
        heavy_file (str): heavy chain input file.
        light_file (str): light chain input file.
        out_file (str): heavy chain output file.
        doublets (str): method for handling multiple heavy chains per cell. one of 'drop' or 'count'.
        format (str): file format. one of 'changeo' or 'airr'.
        """
        # Set column names
        if fileformat == 'changeo':
            cell_id = 'cell_id'
            clone_id = 'clone_id'
            v_call = 'v_call'
            j_call = 'j_call'
            junction_length = 'junction_length'
            umi_count = 'umicount'
        elif fileformat == 'airr':
            cell_id = 'cell_id'
            clone_id = 'clone_id'
            v_call = 'v_call'
            j_call = 'j_call'
            junction_length = 'junction_length'
            umi_count = 'umi_count'
        else:
            sys.exit("Invalid format %s" % fileformat)
    
        # read in heavy and light DFs
        heavy_df = pd.read_csv(heavy_file, dtype='object', na_values=['', 'None', 'NA'], sep='\t')
        light_df = pd.read_csv(light_file, dtype='object', na_values=['', 'None', 'NA'], sep='\t')
    
        # column checking
        expected_heavy_columns = [cell_id, clone_id, v_call, j_call, junction_length, umi_count]
        if set(expected_heavy_columns).issubset(heavy_df.columns) is False:
            raise ValueError("Missing one or more columns in heavy chain file: " + ", ".join(expected_heavy_columns))
        expected_light_columns = [cell_id, v_call, j_call, junction_length, umi_count]
        if set(expected_light_columns).issubset(light_df.columns) is False:
            raise ValueError("Missing one or more columns in light chain file: " + ", ".join(expected_light_columns))
    
        # Fix types
        heavy_df[junction_length] = heavy_df[junction_length].astype('int')
        light_df[junction_length] = light_df[junction_length].astype('int')
    
        # filter multiple heavy chains
        if doublets == 'drop':
            heavy_df = heavy_df.drop_duplicates(cell_id, keep=False)
            if heavy_df.empty is True:
                raise ValueError("Empty heavy chain data, after doublets drop. Are you combining experiments in a single file? If so, split your data into multiple files.")
        elif doublets == 'count':
            heavy_df[umi_count] = heavy_df[umi_count].astype('int')
            heavy_df = heavy_df.groupby(cell_id, sort=False).apply(lambda x: x.nlargest(1, umi_count))
    
        # transfer clone IDs from heavy chain df to light chain df
        clone_dict = {v[cell_id]:v[clone_id] for k, v in heavy_df[[clone_id, cell_id]].T.to_dict().items()}
        light_df = light_df.loc[light_df[cell_id].apply(lambda x: x in clone_dict.keys()), ]
        light_df[clone_id] = light_df.apply(lambda row: clone_dict[row[cell_id]], axis = 1)
    
        # generate a "cluster_dict" of CELL:CLONE dictionary from light df  (TODO: use receptor object V/J gene names)
        cluster_dict = clusterLinkage(light_df[cell_id],
                                    light_df.apply(lambda row:
                                                    getGene(row[v_call]) + ',' + \
                                                    getGene(row[j_call]) + ',' + \
                                                    str(row[junction_length]) + ',' + row[clone_id], axis=1))
    
        # add assignments to heavy_df
        heavy_df = heavy_df.loc[heavy_df[cell_id].apply(lambda x: x in cluster_dict.keys()), :]
        heavy_df[clone_id] = heavy_df[clone_id] + '_' + heavy_df.apply(lambda row: str(cluster_dict[row[cell_id]]), axis=1)
        
        # write heavy chains
        heavy_df.to_csv(out_file, sep='\t', index=False)
        return(heavy_df, light_df)

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd)
    
    h_df, l_df = _lightCluster(h_file2, l_file, outfile, doublets=doublets, fileformat=fileformat)

    h_df = load_data(h_df)
    # create a dictionary for cell_id : clone_id from h_df
    linked_clones = dict(zip(h_df['cell_id'], h_df['clone_id']))
    
    # create a clone_reference 
    clone_ref = list(set(h_df['clone_id']))
    clone_ref = [c.split('_')[1] if c is not np.nan else c for c in clone_ref]
    l_df = load_data(l_df)

    for x in l_df.index:
        if l_df.loc[x, 'clone_id'] in clone_ref:
            l_df.loc[x, 'clone_id'] = linked_clones[l_df.loc[x, 'cell_id']]
        else:
            try:
                l_df.loc[x, 'clone_id'] = l_df.loc[x, 'cell_id']+'_notlinked'
            except:
                pass

    cloned_ = pd.concat([h_df, l_df])
    # transfer the new clone_id to the heavy + light file
    dat['clone_id'] = pd.Series(cloned_['clone_id'])

    if self.__class__ == Dandelion:
        self.__init__(data = dat)
    else:
        out = Dandelion(data = dat)
        return(out)
