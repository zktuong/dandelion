"""Trajectory functions."""
# Created on Mon Sep 19 21:30:44 2022
# @author: chenqu, kp9, kelvin
import re
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp

from collections import Counter
from anndata import AnnData
from typing import List, Optional, Union

from dandelion.utilities._utilities import bh, Literal


def setup_vdj_pseudobulk(
    adata: AnnData,
    mode: Optional[Literal["B", "abT", "gdT"]] = "abT",
    subsetby: Optional[str] = None,
    groups: Optional[List[str]] = None,
    allowed_chain_status: Optional[List[str]] = [
        "Single pair",
        "Extra pair",
        "Extra pair-exception",
        "Orphan VDJ",
        "Orphan VDJ-exception",
    ],
    productive_vdj: bool = True,
    productive_vj: bool = True,
    extract_cols: Optional[List[str]] = None,
    productive_cols: Optional[List[str]] = None,
    check_vdj_mapping: Optional[List[str]] = [
        "v_call",
        "j_call",
    ],
    check_vj_mapping: Optional[List[str]] = [
        "v_call",
        "j_call",
    ],
    check_extract_cols_mapping: Optional[List[str]] = None,
) -> AnnData:
    """Function for prepare anndata for computing pseudobulk vdj feature space.

    Parameters
    ----------
    adata : AnnData
        cell adata before constructing anndata.
    mode : Literal['B', 'abT', 'gdT'], optional
        Optional mode for extractin the V(D)J genes. If set as `None`, it requires the option `extract_cols` to be
        specified with a list of column names where this will be used to retrieve the main call.
    subsetby : Optional[str], optional
        If provided, only the groups/categories in this column will be used for computing the VDJ feature space.
    groups : Optional[List], optional
        If provided, only the following groups/categories will be used for computing the VDJ feature space.
    allowed_chain_status : Optional[List], optional
        If provided, only the ones in this list are kept from the `chain_status` column.
        Defaults to ["Single pair", "Extra pair", "Extra pair-exception", "Orphan VDJ", "Orphan VDJ-exception"].
    productive_vdj: bool, optional
        If True, cells will only be kept if the main VDJ chain is productive.
    productive_vj: bool, optional
        If True, cells will only be kept if the main VJ chain is productive.
    extract_cols : Optional[List[str]], optional
        Column names where VDJ/VJ information is stored so that this will be used instead of the standard columns.
    productive_cols: Optional[List[str]], optional
        Column names where contig productive status is stored so that this will be used instead of the standard columns.
    check_vdj_mapping: Optional[List[str]], optional
        Only columns in the argument will be checked for unclear mapping (containing comma) in VDJ calls.
        Specifying None will skip this step.
    check_vj_mapping: Optional[List[str]], optional
        Only columns in the argument will be checked for unclear mapping (containing comma) in VJ calls.
        Specifying None will skip this step.
    check_vj_extract_cols_mapping: Optional[List[str]], optional
        Only columns in the argument will be checked for unclear mapping (containing comma) in columns specified in extract_cols.
        Specifying None will skip this step.
    Returns
    -------
    AnnData
        filtered cell adata object.
    """
    # keep ony relevant cells (ones with a pair of chains) based on productive column
    if mode is not None:
        if productive_vdj:
            adata = adata[
                adata.obs["productive_" + mode + "_VDJ"].str.startswith("T")
            ].copy()
        if productive_vj:
            adata = adata[
                adata.obs["productive_" + mode + "_VJ"].str.startswith("T")
            ].copy()
    else:
        if productive_cols is not None:
            for col in productive_cols:
                adata = adata[adata.obs[col].str.startswith("T")].copy()

    if allowed_chain_status is not None:
        adata = adata[
            adata.obs["chain_status"].isin(allowed_chain_status)
        ].copy()

    if (groups is not None) and (subsetby is not None):
        adata = adata[adata.obs[subsetby].isin(groups)].copy()

    if "v_call_genotyped_VDJ" in adata.obs:
        v_call = "v_call_genotyped_"
    else:
        v_call = "v_call_"

    if mode is not None:
        adata.obs["v_call_" + mode + "_VDJ_main"] = [
            x.split("|")[0] if x != "None" else "None"
            for x in adata.obs[v_call + mode + "_VDJ"]
        ]
        adata.obs["d_call_" + mode + "_VDJ_main"] = [
            x.split("|")[0] if x != "None" else "None"
            for x in adata.obs["d_call_" + mode + "_VDJ"]
        ]
        adata.obs["j_call_" + mode + "_VDJ_main"] = [
            x.split("|")[0] if x != "None" else "None"
            for x in adata.obs["j_call_" + mode + "_VDJ"]
        ]
        adata.obs["v_call_" + mode + "_VJ_main"] = [
            x.split("|")[0] if x != "None" else "None"
            for x in adata.obs[v_call + mode + "_VJ"]
        ]
        adata.obs["j_call_" + mode + "_VJ_main"] = [
            x.split("|")[0] if x != "None" else "None"
            for x in adata.obs["j_call_" + mode + "_VJ"]
        ]
    else:
        for col in extract_cols:
            adata.obs[col + "_main"] = [
                x.split("|")[0] if x != "None" else "None"
                for x in adata.obs[col]
            ]

    # remove any cells if there's unclear mapping
    if mode is not None:
        if check_vdj_mapping is not None:
            if not isinstance(check_vdj_mapping, list):
                check_vdj_mapping = [check_vdj_mapping]
            for col in check_vdj_mapping:
                adata = adata[
                    ~(
                        adata.obs[col + "_" + mode + "_VDJ_main"].str.contains(
                            ",|None|No_contig"
                        )
                    )
                ].copy()
        if check_vj_mapping is not None:
            if not isinstance(check_vj_mapping, list):
                check_vj_mapping = [check_vj_mapping]
            for col in check_vj_mapping:
                adata = adata[
                    ~(
                        adata.obs[col + "_" + mode + "_VJ_main"].str.contains(
                            ",|None|No_contig"
                        )
                    )
                ].copy()
    else:
        if check_extract_cols_mapping is not None:
            for col in check_extract_cols_mapping:
                adata = adata[
                    ~(adata.obs[col + "_main"].str.contains(",|None|No_contig"))
                ].copy()

    return adata


def _get_pbs(pbs, obs_to_bulk, adata):
    """
    Helper function to ensure we have a cells by pseudobulks matrix which we can use for
    pseudobulking. Uses the pbs and obs_to_bulk inputs to vdj_pseudobulk() and
    gex_pseudobulk().
    """
    # well, we need some way to pseudobulk
    if pbs is None and obs_to_bulk is None:
        raise ValueError(
            "You need to specify `pbs` or `obs_to_bulk` when calling the function"
        )

    # but just one
    if pbs is not None and obs_to_bulk is not None:
        raise ValueError("You need to specify `pbs` or `obs_to_bulk`, not both")

    # turn the pseubodulk matrix dense if need be
    if pbs is not None:
        if sp.sparse.issparse(pbs):
            pbs = pbs.todense()

    # get the obs-derived pseudobulk
    if obs_to_bulk is not None:
        if type(obs_to_bulk) is list:
            # this will create a single value by pasting all the columns together
            tobulk = adata.obs[obs_to_bulk].T.astype(str).agg(",".join)
        else:
            # we just have a single column
            tobulk = adata.obs[obs_to_bulk]
        # this pandas function creates the exact pseudobulk assignment we want
        # this needs to be different than the default uint8
        # as you can have more than 255 cells in a pseudobulk, it turns out
        pbs = pd.get_dummies(tobulk, dtype="uint16").values
    return pbs


def _get_pbs_obs(pbs, obs_to_take, adata):
    """
    Helper function to create the pseudobulk object's obs. Uses the pbs and obs_to_take
    inputs to vdj_pseudobulk() and gex_pseudobulk().
    """
    # prepare per-pseudobulk calls of specified metadata columns
    pbs_obs = pd.DataFrame(index=np.arange(pbs.shape[1]))
    if obs_to_take is not None:
        # just in case a single is passed as a string
        if type(obs_to_take) is not list:
            obs_to_take = [obs_to_take]
        # now we can iterate over this nicely
        # using the logic of milopy's annotate_nhoods()
        for anno_col in obs_to_take:
            anno_dummies = pd.get_dummies(adata.obs[anno_col])
            # this needs to be turned to a matrix so dimensions get broadcast correctly
            anno_count = np.asmatrix(pbs).T.dot(anno_dummies.values)
            anno_frac = np.array(anno_count / anno_count.sum(1))
            anno_frac = pd.DataFrame(
                anno_frac,
                index=np.arange(pbs.shape[1]),
                columns=anno_dummies.columns,
            )
            pbs_obs[anno_col] = anno_frac.idxmax(1)
            pbs_obs[anno_col + "_fraction"] = anno_frac.max(1)
    # report the number of cells for each pseudobulk
    # ensure pbs is an array so that it sums into a vector that can go in easily
    pbs_obs["cell_count"] = np.sum(np.asarray(pbs), axis=0)
    return pbs_obs


def vdj_pseudobulk(
    adata: AnnData,
    pbs: Optional[Union[np.ndarray, sp.sparse.csr_matrix]] = None,
    obs_to_bulk: Optional[Union[str, List[str]]] = None,
    obs_to_take: Optional[Union[str, List[str]]] = None,
    cols: Optional[List[str]] = None,
) -> AnnData:
    """Function for making pseudobulk vdj feature space. One of `pbs` or `obs_to_bulk`
    needs to be specified when calling.

    Parameters
    ----------
    adata : AnnData
        Cell adata, preferably after `ddl.tl.setup_vdj_pseudobulk()`
    pbs: Optional[array], optional
        Optional binary matrix with cells as rows and pseudobulk groups as columns
    obs_to_bulk: Optional[Union[str, List[str]]], optional
        Optional obs column(s) to group pseudobulks into; if multiple are provided, they
        will be combined
    obs_to_take: Optional[Union[str, List[str]]]
        Optional obs column(s) to identify the most common value of for each pseudobulk.
    cols: Optional[List], optional
        If provided, use the specified obs columns to extract V(D)J calls

    Returns
    -------
    AnnData
        pb_adata, whereby each observation is a pseudobulk:\n
        VDJ usage frequency stored in pb_adata.X\n
        VDJ genes stored in pb_adata.var\n
        pseudobulk metadata stored in pb_adata.obs\n
        pseudobulk assignment (binary matrix with input cells as columns) stored in pb_adata.obsm['pbs']\n
    """
    # get our cells by pseudobulks matrix
    pbs = _get_pbs(pbs, obs_to_bulk, adata)

    # if not specified by the user, use the following default dandelion VJ columns
    if cols is None:
        cols = [i for i in adata.obs if re.search("_VDJ_main|_VJ_main", i)]

    # perform matrix multiplication of pseudobulks by cells matrix by a cells by VJs matrix
    # start off by creating the cell by VJs matrix, using the pandas dummies again
    # skip the prefix stuff as the VJ genes will be unique in the columns
    vjs = pd.get_dummies(adata.obs[cols], prefix="", prefix_sep="")
    # TODO: DENAN SOMEHOW? AS IN NAN GENES?
    # can now multiply transposed pseudobulk assignments by this vjs thing, turn to df
    vj_pb_count = pbs.T.dot(vjs.values)
    df = pd.DataFrame(
        vj_pb_count, index=np.arange(pbs.shape[1]), columns=vjs.columns
    )

    # loop over V(D)J gene categories
    for col in cols:
        # identify columns holding genes belonging to the category
        # and then normalise the values to 1 for each pseudobulk
        mask = np.isin(df.columns, np.unique(adata.obs[col]))
        df.loc[:, mask] = df.loc[:, mask].div(
            df.loc[:, mask].sum(axis=1), axis=0
        )

    # create obs for the new pseudobulk object
    pbs_obs = _get_pbs_obs(pbs, obs_to_take, adata)

    # store our feature space and derived metadata into an AnnData
    pb_adata = sc.AnnData(
        np.array(df), var=pd.DataFrame(index=df.columns), obs=pbs_obs
    )
    # store the pseudobulk assignments, as a sparse for storage efficiency
    # transpose as the original matrix is cells x pseudobulks
    pb_adata.obsm["pbs"] = sp.sparse.csr_matrix(pbs.T)
    return pb_adata


def pseudotime_transfer(
    adata: AnnData, pr_res: "palantir.presults.PResults", suffix: str = ""
):
    """Function to add pseudotime and branch probabilities into adata.obs in place.

    Parameters
    ----------
    adata : AnnData
        adata for which pseudotime to be transferred to
    pr_res : palantir.presults.PResults
        palantir pseudotime inference output object
    suffix : str, optional
        suffix to be added after the added column names, default "" (none)
    """
    adata.obs["pseudotime" + suffix] = pr_res.pseudotime.copy()

    for col in pr_res.branch_probs.columns:
        adata.obs["prob_" + col + suffix] = pr_res.branch_probs[col].copy()

    return adata


def project_pseudotime_to_cell(
    adata: AnnData, pb_adata: AnnData, term_states: List[str], suffix: str = ""
) -> AnnData:
    """Function to project pseudotime & branch probabilities from pb_adata (pseudobulk adata) to adata (cell adata).

    Parameters
    ----------
    adata : AnnData
        Cell adata, preferably after `ddl.tl.setup_vdj_pseudobulk()`
    pb_adata : AnnData
        neighbourhood/pseudobulked adata
    term_states : List[str]
        list of terminal states with branch probabilities to be transferred
    suffix : str, optional
        suffix to be added after the added column names, default "" (none)

    Returns
    -------
    AnnData
        subset of adata whereby cells that don't belong to any neighbourhood are removed
        and projected pseudotime information stored in .obs - `pseudotime+suffix`, and `'prob_'+term_state+suffix` for each terminal state
    """
    # extract out cell x pseudobulk matrix. it's stored as pseudobulk x cell so transpose
    nhoods = np.array(pb_adata.obsm["pbs"].T.todense())

    # leave out cells that don't belong to any neighbourhood
    nhoodsum = np.sum(nhoods, axis=1)
    cdata = adata[nhoodsum > 0].copy()
    print(
        "number of cells removed due to not belonging to any neighbourhood",
        sum(nhoodsum == 0),
    )  # print how many cells removed
    # also subset the pseudobulk_assignments
    pb_assign_trim = nhoods[nhoodsum > 0]

    # for each cell pseudotime_mean is the average of the pseudotime of all pseudobulks the cell is in, weighted by 1/neighbourhood size
    nhoods_cdata = nhoods[nhoodsum > 0, :]
    nhoods_cdata_norm = nhoods_cdata / np.sum(
        nhoods_cdata, axis=0, keepdims=True
    )

    col_list = ["pseudotime" + suffix] + [
        "prob_" + state + suffix for state in term_states
    ]
    for col in col_list:
        cdata.obs[col] = (
            np.array(
                nhoods_cdata_norm.dot(pb_adata.obs[col]).T
                / np.sum(nhoods_cdata_norm, axis=1)
            )
            .flatten()
            .copy()
        )
    cdata.uns["pseudobulk_assignments"] = pb_assign_trim.copy()

    return cdata


def pseudobulk_gex(
    adata_raw: AnnData,
    pbs: Optional[Union[np.ndarray, sp.sparse.csr_matrix]] = None,
    obs_to_bulk: Optional[Union[str, List[str]]] = None,
    obs_to_take: Optional[Union[str, List[str]]] = None,
) -> AnnData:
    """Function to pseudobulk gene expression (raw count).

    Parameters
    ----------
    adata_raw : AnnData
        Needs to have raw counts in .X
    pbs: Optional[array], optional
        Optional binary matrix with cells as rows and pseudobulk groups as columns
    obs_to_bulk: Optional[Union[str, List[str]]], optional
        Optional obs column(s) to group pseudobulks into; if multiple are provided, they
        will be combined
    obs_to_take: Optional[Union[str, List[str]]]
        Optional obs column(s) to identify the most common value of for each pseudobulk

    Returns
    -------
    AnnData
        pb_adata whereby each observation is a cell neighbourhood\n
        pseudobulked gene expression stored in pb_adata.X\n
        genes stored in pb_adata.var\n
        pseudobulk metadata stored in pb_adata.obs\n
        pseudobulk assignment (binary matrix with input cells as columns) stored in pb_adata.obsm['pbs']\n
    """
    # get our cells by pseudobulks matrix
    pbs = _get_pbs(pbs, obs_to_bulk, adata_raw)

    # make pseudobulk matrix
    pbs_X = adata_raw.X.T.dot(pbs)

    # create obs for the new pseudobulk object
    pbs_obs = _get_pbs_obs(pbs, obs_to_take, adata_raw)

    ## Make new anndata object
    pb_adata = sc.AnnData(pbs_X.T, obs=pbs_obs, var=adata_raw.var)
    # store the pseudobulk assignments, as a sparse for storage efficiency
    # transpose as the original matrix is cells x pseudobulks
    pb_adata.obsm["pbs"] = sp.sparse.csr_matrix(pbs.T)
    return pb_adata


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
