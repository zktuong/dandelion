import re
import scipy.stats

import numpy as np
import pandas as pd
import scanpy as sc

from anndata import AnnData
from sklearn import mixture

from dandelion.utilities._utilities import bh


def recipe_scanpy_qc(
    adata: AnnData,
    layer: str | None = None,
    mito_startswith: str = "MT-",
    max_genes: int = 2500,
    min_genes: int = 200,
    mito_cutoff: int | None = 5,
    run_scrublet: bool = True,
    pval_cutoff: float = 0.1,
    min_counts: int | None = None,
    max_counts: int | None = None,
    blacklist: list[str] | None = None,
    vdj_pattern: str = "^TR[AB][VDJ]|^IG[HKL][VDJC]",
):
    """
    Recipe for running a standard scanpy QC workflow.

    Parameters
    ----------
    adata : AnnData
        annotated data matrix of shape n_obs × n_vars. Rows correspond to cells
        and columns to genes.
    layer : str | None, optional
        name of layer to run scrublet on if supplied.
    mito_startswith : str, optional
        string pattern used for searching mitochondrial genes.
    max_genes : int, optional
        maximum number of genes expressed required for a cell to pass filtering
    min_genes : int, optional
        minimum number of genes expressed required for a cell to pass filtering
    mito_cutoff : int | None, optional
        maximum percentage mitochondrial content allowed for a cell to pass filtering.
    run_scrublet : bool, optional
        whether or not to run scrublet for doublet detection.
    pval_cutoff : float, optional
        maximum Benjamini-Hochberg corrected p value from doublet detection
        protocol allowed for a cell to pass filtering. Default is 0.1.
    min_counts : int | None, optional
        minimum number of counts required for a cell to pass filtering.
    max_counts : int | None, optional
        maximum number of counts required for a cell to pass filtering.
    blacklist : list[str] | None, optional
        if provided, will exclude these genes from highly variable genes list.
    vdj_pattern : str, optional
        string pattern for search VDJ genes to exclude from highly variable genes.

    Raises
    ------
    ImportError
        if `scrublet` not installed.
    """
    _adata = adata.copy()
    # run basic scanpy pipeline
    sc.pp.filter_cells(_adata, min_genes=0)
    _adata.var["mt"] = _adata.var_names.str.startswith(mito_startswith)
    sc.pp.calculate_qc_metrics(
        _adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    if mito_cutoff is None:
        # use a model-based method to determine the cut off
        # for mitochondrial content
        gmm = mixture.GaussianMixture(
            n_components=2, max_iter=1000, covariance_type="full"
        )
        X = _adata.obs[["pct_counts_mt", "n_genes_by_counts"]]
        _adata.obs["gmm_pct_count_clusters"] = gmm.fit(X).predict(X)
        # use a simple metric to workout which cluster
        # is the one that contains lower mito content?
        A1 = (
            _adata[_adata.obs["gmm_pct_count_clusters"] == 0]
            .obs["pct_counts_mt"]
            .mean()
        )
        B1 = (
            _adata[_adata.obs["gmm_pct_count_clusters"] == 1]
            .obs["pct_counts_mt"]
            .mean()
        )
        A2 = (
            _adata[_adata.obs["gmm_pct_count_clusters"] == 0]
            .obs["n_genes_by_counts"]
            .mean()
        )
        B2 = (
            _adata[_adata.obs["gmm_pct_count_clusters"] == 1]
            .obs["n_genes_by_counts"]
            .mean()
        )
        if (A1 > B1) and (A2 < B2):
            keepdict = {0: False, 1: True}
        else:
            keepdict = {1: False, 0: True}
        _adata.obs["gmm_pct_count_clusters_keep"] = [
            keepdict[x] for x in _adata.obs["gmm_pct_count_clusters"]
        ]
    if run_scrublet:
        # run scrublet
        try:
            import scrublet as scr
        except ImportError:
            raise ImportError(
                "Please install scrublet with pip install scrublet."
            )
        if layer is None:
            scrub = scr.Scrublet(_adata.X)
        else:
            scrub = scr.Scrublet(_adata.layers[layer])
        doublet_scores, _ = scrub.scrub_doublets(verbose=False)
        _adata.obs["scrublet_score"] = doublet_scores
        # overcluster prep.
        sc.pp.normalize_total(_adata, target_sum=1e4)
        sc.pp.log1p(_adata)
        sc.pp.highly_variable_genes(
            _adata, min_mean=0.0125, max_mean=3, min_disp=0.5
        )
        for i in _adata.var.index:
            if vdj_pattern is not None:
                if re.search(vdj_pattern, i):
                    _adata.var.at[i, "highly_variable"] = False
            if blacklist is not None:
                if i in blacklist:
                    _adata.var.at[i, "highly_variable"] = False
        _adata = _adata[:, _adata.var["highly_variable"]].copy()
        sc.pp.scale(_adata, max_value=10)
        sc.tl.pca(_adata, svd_solver="arpack")
        sc.pp.neighbors(_adata, n_neighbors=10, n_pcs=50)
        # overclustering proper - do basic clustering first,
        # then cluster each cluster
        sc.tl.leiden(_adata)
        for clus in list(np.unique(_adata.obs["leiden"]))[0]:
            sc.tl.leiden(
                _adata, restrict_to=("leiden", [clus]), key_added="leiden_R"
            )
        # weird how the new anndata/scanpy is forcing this
        for clus in list(np.unique(_adata.obs["leiden"]))[1:]:
            sc.tl.leiden(
                _adata, restrict_to=("leiden_R", [clus]), key_added="leiden_R"
            )
        # compute the cluster scores - the median of Scrublet scores per
        # overclustered cluster
        for clus in np.unique(_adata.obs["leiden_R"]):
            _adata.obs.loc[
                _adata.obs["leiden_R"] == clus, "scrublet_cluster_score"
            ] = np.median(
                _adata.obs.loc[_adata.obs["leiden_R"] == clus, "scrublet_score"]
            )
        # now compute doublet p-values. figure out the median and mad
        # (from above-median values) for the distribution
        med = np.median(_adata.obs["scrublet_cluster_score"])
        mask = _adata.obs["scrublet_cluster_score"] > med
        mad = np.median(_adata.obs["scrublet_cluster_score"][mask] - med)
        # 1 sided test for catching outliers
        pvals = 1 - scipy.stats.norm.cdf(
            _adata.obs["scrublet_cluster_score"], loc=med, scale=1.4826 * mad
        )
        _adata.obs["scrublet_score_bh_pval"] = bh(pvals)
        # threshold the p-values to get doublet calls.
        _adata.obs["is_doublet"] = (
            _adata.obs["scrublet_score_bh_pval"] < pval_cutoff
        )
    else:
        _adata.obs["is_doublet"] = False
    if mito_cutoff is not None:
        if min_counts is None and max_counts is None:
            _adata.obs["filter_rna"] = (
                (
                    pd.Series(
                        [
                            ((n < min_genes) or (n > max_genes))
                            for n in _adata.obs["n_genes_by_counts"]
                        ],
                        index=_adata.obs.index,
                    )
                )
                | (_adata.obs["pct_counts_mt"] >= mito_cutoff)
                | (_adata.obs.is_doublet)
            )
        else:
            if min_counts is not None:
                if max_counts is not None:
                    _adata.obs["filter_rna"] = (
                        (
                            pd.Series(
                                [
                                    ((n < min_genes) or (n > max_genes))
                                    for n in _adata.obs["n_genes_by_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | (
                            pd.Series(
                                [
                                    min_counts < n > max_counts
                                    for n in _adata.obs["total_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | (_adata.obs["pct_counts_mt"] >= mito_cutoff)
                        | (_adata.obs.is_doublet)
                    )
                else:
                    _adata.obs["filter_rna"] = (
                        (
                            pd.Series(
                                [
                                    ((n < min_genes) or (n > max_genes))
                                    for n in _adata.obs["n_genes_by_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | (
                            pd.Series(
                                [
                                    n < min_counts
                                    for n in _adata.obs["total_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | (_adata.obs["pct_counts_mt"] >= mito_cutoff)
                        | (_adata.obs.is_doublet)
                    )
            else:
                if max_counts is not None:
                    _adata.obs["filter_rna"] = (
                        (
                            pd.Series(
                                [
                                    ((n < min_genes) or (n > max_genes))
                                    for n in _adata.obs["n_genes_by_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | (
                            pd.Series(
                                [
                                    n > max_counts
                                    for n in _adata.obs["total_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | (_adata.obs["pct_counts_mt"] >= mito_cutoff)
                        | (_adata.obs.is_doublet)
                    )
    else:
        if min_counts is None and max_counts is None:
            _adata.obs["filter_rna"] = (
                (
                    pd.Series(
                        [
                            ((n < min_genes) or (n > max_genes))
                            for n in _adata.obs["n_genes_by_counts"]
                        ],
                        index=_adata.obs.index,
                    )
                )
                | ~(_adata.obs.gmm_pct_count_clusters_keep)
                | (_adata.obs.is_doublet)
            )
        else:
            if min_counts is not None:
                if max_counts is not None:
                    _adata.obs["filter_rna"] = (
                        (
                            pd.Series(
                                [
                                    ((n < min_genes) or (n > max_genes))
                                    for n in _adata.obs["n_genes_by_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | (
                            pd.Series(
                                [
                                    min_counts < n > max_counts
                                    for n in _adata.obs["total_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | ~(_adata.obs.gmm_pct_count_clusters_keep)
                        | (_adata.obs.is_doublet)
                    )
                else:
                    _adata.obs["filter_rna"] = (
                        (
                            pd.Series(
                                [
                                    ((n < min_genes) or (n > max_genes))
                                    for n in _adata.obs["n_genes_by_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | (
                            pd.Series(
                                [
                                    n < min_counts
                                    for n in _adata.obs["total_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | ~(_adata.obs.gmm_pct_count_clusters_keep)
                        | (_adata.obs.is_doublet)
                    )
            else:
                if max_counts is not None:
                    _adata.obs["filter_rna"] = (
                        (
                            pd.Series(
                                [
                                    ((n < min_genes) or (n > max_genes))
                                    for n in _adata.obs["n_genes_by_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | (
                            pd.Series(
                                [
                                    n > max_counts
                                    for n in _adata.obs["total_counts"]
                                ],
                                index=_adata.obs.index,
                            )
                        )
                        | ~(_adata.obs.gmm_pct_count_clusters_keep)
                        | (_adata.obs.is_doublet)
                    )
    bool_dict = {True: "True", False: "False"}

    _adata.obs["is_doublet"] = [bool_dict[x] for x in _adata.obs["is_doublet"]]
    _adata.obs["filter_rna"] = [bool_dict[x] for x in _adata.obs["filter_rna"]]

    # removing columns that probably don't need anymore
    if mito_cutoff is not None:
        _adata.obs = _adata.obs.drop(
            [
                "leiden",
                "leiden_R",
                "scrublet_cluster_score",
                "scrublet_score_bh_pval",
            ],
            axis=1,
        )
    else:
        _adata.obs = _adata.obs.drop(
            [
                "leiden",
                "leiden_R",
                "scrublet_cluster_score",
                "scrublet_score_bh_pval",
                "gmm_pct_count_clusters",
            ],
            axis=1,
        )
    if not run_scrublet:
        _adata.obs = _adata.obs.drop(["is_doublet"], axis=1)
    adata.obs = _adata.obs.copy()
