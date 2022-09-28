"""Trajectory functions."""
# Created on Mon Sep 19 21:30:44 2022
# @author: chenqu

import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
import palantir

from collections import Counter
from anndata import AnnData
from typing import List

from ..utilities._utilities import bh


def vdj_nhood(
    adata: AnnData,
    cols: List[str] = [
        "ab_IR_VDJ_1_v_gene",
        "ab_IR_VDJ_1_j_gene",
        "ab_IR_VJ_1_v_gene",
        "ab_IR_VJ_1_j_gene",
    ],
) -> AnnData:
    """Function for making neighbourhood vdj feature space.

    Parameters
    ----------
    adata : AnnData
        cell adata with neighbourhood information stored in adata.uns['nhood_adata'] & adata.obsm['nhoods']
    cols : List[str], optional
        columns of VDJ in adata.obs to be used here i.e. columns of e.g. TRAV, TRAJ, TRBV, TRBJ

    Returns
    -------
    AnnData
        nhood_adata, whereby each observation is a cell neighbourhood:
        VDJ usage frequency stored in nhood_adata.X
        VDJ genes stored in nhood_adata.var
        neighbourhood metadata stored in nhood_adata.obs
    """
    nhoods = adata.obsm["nhoods"].todense()
    # the .ravel() turns the four V(D)J columns to a single vector of values
    # and doing a .unique() of that gets us all the possible genes present in the object
    # we want a column for every single V(D)J gene encountered, so this is perfect
    df = pd.DataFrame(
        0,
        columns=pd.unique(adata.obs[cols].values.ravel("K")),
        index=np.arange(nhoods.shape[1]),
    )
    for i in df.index:
        # extract the metadata for this neighbourhood, and just the V(D)J gene columns
        sub = adata.obs.loc[nhoods[:, i] == 1, cols]
        # quickly count up the .ravel()ed form up, and this goes straight into a df - neat!
        df.loc[i, :] = Counter(sub.values.ravel("K"))
    # the above procedure leaves NaNs in unencountered locations
    df.fillna(0, inplace=True)
    # TODO: DENAN SOMEHOW? AS IN NAN GENES?
    # loop over V(D)J gene categories
    for col in cols:
        # identify columns holding genes belonging to the category
        # and then normalise the values to 1 for each neighbourhood
        mask = np.isin(df.columns, np.unique(adata.obs[col]))
        df.loc[:, mask] = df.loc[:, mask].div(
            df.loc[:, mask].sum(axis=1), axis=0
        )
    # store our feature space and some existing metadata into an AnnData
    nhood_adata = sc.AnnData(
        np.array(df),
        var=pd.DataFrame(index=df.columns),
        obs=adata.uns["nhood_adata"].obs,
    )
    return nhood_adata


def pseudotime_transfer(
    adata: AnnData, pr_res: palantir.presults.PResults, suffix: str
):
    """Function to add pseudotime and branch probabilities into adata.obs in place.

    Parameters
    ----------
    adata : AnnData
        nhood_adata for which pseudotime to be transferred to
    pr_res : palantir.presults.PResults
        palantir pseudotime inference output object
    suffix : str
        suffix to be added after the added column names
    """
    adata.obs["pseudotime" + suffix] = pr_res.pseudotime.copy()

    for col in pr_res.branch_probs.columns:
        adata.obs["prob_" + col + suffix] = pr_res.branch_probs[col].copy()


def pseudotime_cell(
    adata: AnnData, nhood_adata: AnnData, term_states: List[str], suffix: str
) -> AnnData:
    """Function to project pseudotime & branch probabilities from nhood_adata (neighbourhood adata) to adata (cell adata).

    Parameters
    ----------
    adata : AnnData
        cell adata
    nhood_adata : AnnData
        neighbourhood adata
    term_states : List[str]
        list of terminal states with branch probabilities to be transferred
    suffix : str
        suffix to be added after the added column names

    No Longer Returned
    ------------------
    TYPE
        subset of adata whereby cells that don't belong to any neighbourhood are removed
        and projected pseudotime information stored in .obs - ['pseudotime'+suffix], and ['prob_'+term_state+suffix] for each terminal state
    """
    # extract out cell x neighbourhood matrix
    nhoods = np.array(adata.obsm["nhoods"].todense())

    # leave out cells that don't belong to any neighbourhood
    cdata = adata[np.sum(nhoods, axis=1) > 0].copy()
    print(
        "number of cells removed", sum(np.sum(nhoods, axis=1) == 0)
    )  # print how many cells removed

    # for each cell pseudotime_mean is the average of the pseudotime of all neighbourhoods the cell is in, weighted by 1/neighbourhood size
    nhoods_cdata = nhoods[np.sum(nhoods, axis=1) > 0, :]
    nhoods_cdata_norm = nhoods_cdata / np.sum(
        nhoods_cdata, axis=0, keepdims=True
    )

    col_list = ["pseudotime" + suffix] + [
        "prob_" + state + suffix for state in term_states
    ]
    for col in col_list:
        cdata.obs[col] = (
            np.array(
                nhoods_cdata_norm.dot(nhood_adata.obs[col]).T
                / np.sum(nhoods_cdata_norm, axis=1)
            )
            .flatten()
            .copy()
        )

    return cdata


def nhood_gex(
    adata: AnnData, adata_raw: AnnData, normalize_log: bool
) -> AnnData:
    """Function to pseudobulk gene expression (raw count) by cell neighbourhoods.

    Parameters
    ----------
    adata : AnnData
        cell adata
    adata_raw : AnnData
        same cells with raw counts. Not normalised, not log1p, just raw counts.
    normalize_log : bool
        True or False, pseudobulked expression to be normalised and log transformed

    Returns
    -------
    AnnData
        nhood_adata whereby each observation is a cell neighbourhood
        pseudobulked gene expression stored in nhood_adata.X
        genes stored in nhood_adata.var
        neighbourhood metadata stored in nhood_adata.obs
    """
    sample_dummies = adata.obsm["nhoods"]

    # make pseudobulk matrix
    pseudobulk_X = adata_raw.X.T.dot(sample_dummies)

    ## Make new anndata object
    nhood_adata = sc.AnnData(
        pseudobulk_X.T, obs=adata.uns["nhood_adata"].obs, var=adata_raw.var
    )

    nhood_adata.raw = nhood_adata.copy()

    # normalise and log1p
    if normalize_log:
        sc.pp.normalize_per_cell(nhood_adata, counts_per_cell_after=10e4)
        sc.pp.log1p(nhood_adata)

    return nhood_adata


def bin_expression(
    adata: AnnData, bin_no: int, genes: List[str], pseudotime_col: str
) -> pd.DataFrame:
    """Function to compute average gene expression in bins along pseudotime.

    Parameters
    ----------
    adata : AnnData
        cell adata.
    bin_no : int
        number of bins to be divided along pseudotime.
    genes : List
        list of genes for the computation
    pseudotime_col : str
        column in adata.obs where pseudotime is stored

    Returns
    -------
    pd.DataFrame
        a dataframe with genes as rows, and pseudotime bins as columns, and averaged gene expression as the data
    """
    # define bins
    bins = np.linspace(0, 1, bin_no + 1)

    # get gene expression
    x = np.array(adata[:, genes].X.todense())
    # get pseudotime
    y = np.array(adata.obs[pseudotime_col])

    # calculate average gene expression in each bin
    gene_summary = pd.DataFrame(columns=bins[:-1], index=genes)
    for i in range(gene_summary.shape[1]):
        time = bins[i]
        select = np.array(bins[i] <= y) & np.array(y < bins[i + 1])
        gene_summary.loc[:, time] = np.mean(x[select, :], axis=0)

    return gene_summary


def chatterjee_corr(
    adata: AnnData, genes: List[str], pseudotime_col: str
) -> pd.DataFrame:
    """Function to compute chatterjee correlation of gene expression with pseudotime.

    Parameters
    ----------
    adata : AnnData
        cell adata
    genes : List[str]
        List of genes selected to compute the correlation
    pseudotime_col : str
        olumn in adata.obs where pseudotime is stored

    Returns
    -------
    pd.DataFrame
        a dataframe with genes as rows, with cor_res (correlation statistics),
        pval (p-value),  adj_pval (p-value adjusted by BH method) as columns.
    """
    # get gene expression
    x = np.array(adata[:, genes].X.todense())
    # add small perturbation for random tie breaking
    x = x + np.random.randn(x.shape[0], x.shape[1]) * 1e-15
    # get pseudotime
    y = list(adata.obs[pseudotime_col])

    # compute chatterjee correlation
    # ref: Sourav Chatterjee (2021) A New Coefficient of Correlation, Journal of the American Statistical Association, 116:536, 2009-2022, DOI: 10.1080/01621459.2020.1758115
    stat = 1 - np.sum(
        np.abs(np.diff(np.argsort(x[np.argsort(y), :], axis=0), axis=0)), axis=0
    ) * 3 / (x.shape[0] ** 2 - 1)
    stat = np.array(stat).flatten()

    pval = 1 - sp.stats.norm.cdf(stat, loc=0, scale=np.sqrt(2 / 5 / x.shape[0]))

    # put results into dataframe cor_res
    cor_res = pd.DataFrame({"cor_stat": stat, "pval": pval})
    cor_res.index = genes

    # compute adjusted pval using BH method
    cor_res["adj_pval"] = bh(cor_res["pval"].to_numpy())

    # sort genes based on adjusted pval
    cor_res = cor_res.sort_values(by="adj_pval")

    return cor_res
