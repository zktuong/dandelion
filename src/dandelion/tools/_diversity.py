#!/usr/bin/env python
import numpy as np
import networkx as nx
import pandas as pd

from anndata import AnnData
from time import sleep
from tqdm import tqdm
from scanpy import logging as logg
from scipy.special import gammaln
from typing import Literal

import dandelion
from dandelion.external.skbio._chao1 import chao1
from dandelion.external.skbio._gini import gini_index
from dandelion.external.skbio._shannon import shannon
from dandelion.tools._network import (
    clone_centrality,
    clone_degree,
    generate_network,
)
from dandelion.utilities._core import *
from dandelion.utilities._io import *
from dandelion.utilities._utilities import *


def clone_rarefaction(
    vdj_data: Dandelion | AnnData,
    groupby: str,
    clone_key: str | None = None,
    rarefaction_key: str | None = None,
    verbose: bool = False,
) -> dict:
    """
    Return rarefaction predictions for cell numbers vs clone size.

    Parameters
    ----------
    vdj_data : Dandelion | AnnData
        Dandelion or AnnData object.
    groupby : str
        Column name to split the calculation of clone numbers for a given number of cells for e.g. sample id, patient etc.
    clone_key : str | None, optional
        Column name specifying the clone_id column in metadata/obs.
    rarefaction_key : str | None, optional
        key for 'diversity' results in AnnData's `.uns`.
    verbose : bool, optional
        whether to print progress.

    Returns
    -------
    dict
        Rarefaction predictions for cell numbers vs clone size.
    """
    if isinstance(vdj_data, AnnData):
        _metadata = vdj_data.obs.copy()
    elif isinstance(vdj_data, Dandelion):
        _metadata = vdj_data.metadata.copy()

    clonekey = clone_key if clone_key is not None else "clone_id"

    groups = list(set(_metadata[groupby]))
    if type(_metadata[clonekey]) == "category":
        _metadata[clonekey] = _metadata[clonekey].cat.remove_unused_categories()
    res = {}
    for g in groups:
        __metadata = _metadata[_metadata[groupby] == g]
        res[g] = __metadata[clonekey].value_counts()
    res_ = pd.DataFrame.from_dict(res, orient="index")

    # remove those with no counts
    logg.info(
        "removing due to zero counts: "
        ", ".join(
            [res_.index[i] for i, x in enumerate(res_.sum(axis=1) == 0) if x]
        ),
    )
    res_ = res_[~(res_.sum(axis=1) == 0)]

    # set up for calculating rarefaction
    tot = res_.apply(sum, axis=1)
    # S = res_.apply(lambda x: x[x > 0].shape[0], axis=1)
    nr = res_.shape[0]

    # append the results to a dictionary
    rarecurve = {}
    sleep(0.5)
    for i in tqdm(
        range(0, nr),
        desc="Calculating rarefaction curve ",
        bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
        disable=not verbose,
    ):
        n = np.arange(1, tot[i], step=10)
        if n[-1:] != tot[i]:
            n = np.append(n, tot[i])
        rarecurve[res_.index[i]] = [
            rarefun(
                np.array(res_.iloc[i,]),
                z,
            )
            for z in n
        ]
    y = pd.DataFrame([rarecurve[c] for c in rarecurve]).T
    pred = pd.DataFrame(
        [np.append(np.arange(1, s, 10), s) for s in res_.sum(axis=1)],
        index=res_.index,
    ).T

    rarefaction_key = (
        rarefaction_key if rarefaction_key is not None else "rarefaction"
    )

    if isinstance(vdj_data, AnnData):
        if rarefaction_key not in vdj_data.uns:
            vdj_data.uns[rarefaction_key] = {}
        vdj_data.uns[rarefaction_key] = {
            "rarefaction_cells_x": pred,
            "rarefaction_clones_y": y,
        }
    if isinstance(vdj_data, Dandelion):
        return {"rarefaction_cells_x": pred, "rarefaction_clones_y": y}


def clone_diversity(
    vdj_data: Dandelion | AnnData,
    groupby: str,
    method: Literal["gini", "chao1", "shannon"] = "gini",
    metric: Literal["clone_network", "clone_degree", "clone_centrality"] = None,
    clone_key: str | None = None,
    min_size: int | None = None,
    n_boot: int = 1000,
    normalize: bool = True,
    use_network: bool = True,
    expanded_only: bool = False,
    use_contracted: bool = False,
    verbose: bool = False,
    **kwargs,
) -> tuple[pd.DataFrame, dict[list[float]]]:
    """
    Compute clonal diversity with bootstrapping.

    Parameters
    ----------
    vdj_data : Dandelion | AnnData
        Dandelion or AnnData object.
    groupby : str
        Column name to calculate the gini indices on, for e.g. sample id, patient etc.
    method : Literal["gini", "chao1", "shannon"], optional
        Method for diversity estimation. Either one of ['gini', 'chao1', 'shannon'].
    metric : Literal["clone_network", "clone_degree", "clone_centrality"], optional
        Metric to use for calculating Gini indices of clones.
        Accepts one of ['clone_network', 'clone_degree', 'clone_centrality'].
        `None` defaults to 'clone_network'.
    clone_key : str | None, optional
        Column name specifying the clone_id column in metadata.
    min_size : int | None, optional
        Minimum cell numbers to keep for diversity calculation. If None, defaults to size of smallest sample.
        Beware that this may lead to very small sample sizes and unreliable estimates if left as None.
    n_boot : int, optional
        Number of times to perform resampling. Default is 50.
    normalize : bool, optional
        Whether or not to return normalized Shannon Entropy according to https://math.stackexchange.com/a/945172.
        Default is True.
    use_network : bool, optional
        Whether or not to use network-based Gini index calculation. Default is True.
    expanded_only : bool, optional
        Whether or not to calculate gini indices using expanded clones only. Default is False i.e. use all cells/clones.
    use_contracted : bool, optional
        Whether or not to perform the gini calculation after contraction of clone network.
        Only applies to calculation of clone size gini index. Default is False.
        This is to try and preserve the single-cell properties of the network.
    verbose : bool, optional
        whether to print progress.
    **kwargs
        passed to dandelion.tl.generate_network

    Returns
    -------
    tuple[pd.DataFrame, dict[list[float]]]
        pandas DataFrame holding summarised diversity estimation and the raw bootstrap results.
    """
    # if AnnData, cannot use network-based gini
    use_network = False if isinstance(vdj_data, AnnData) else use_network
    if (method == "gini") and use_network:
        return diversity_gini(
            vdj_data,
            groupby=groupby,
            metric=metric,
            clone_key=clone_key,
            min_size=min_size,
            n_boot=n_boot,
            expanded_only=expanded_only,
            use_contracted=use_contracted,
            verbose=verbose,
            **kwargs,
        )
    else:
        return diversity_estimates(
            vdj_data,
            method=method,
            groupby=groupby,
            clone_key=clone_key,
            normalize=normalize,
            min_size=min_size,
            n_boot=n_boot,
            verbose=verbose,
        )


def clone_networkstats(
    vdj_data: Dandelion,
    expanded_only: bool = False,
    network_clustersize: bool = False,
    verbose: bool = False,
) -> tuple[defaultdict, defaultdict, defaultdict]:
    """Retrieve network stats.

    Parameters
    ----------
    vdj_data : Dandelion
        input object
    expanded_only : bool, optional
        whether or not to calculate only on expanded clones.
    network_clustersize : bool, optional
        depends on metric.
    verbose : bool, optional
        whether to print progress.

    Returns
    -------
    tuple[defaultdict, defaultdict, defaultdict]
        output nodes names, vertex sizes and cluster sizes.

    Raises
    ------
    AttributeError
        if graph not found.
    TypeError
        if input object is not Dandelion.
    """
    logg.info("Calculating vertex size of nodes after contraction")

    if isinstance(vdj_data, Dandelion):
        if vdj_data.graph is None:
            raise AttributeError(
                "Graph not found. Please run tl.generate_network."
            )
        else:
            if expanded_only:
                G = vdj_data.graph[1]
            else:
                G = vdj_data.graph[0]
            remove_edges = defaultdict(list)
            vertexsizes = defaultdict(list)
            clustersizes = defaultdict(list)
            nodes_names = defaultdict(list)

            for subg in tqdm(
                nx.connected_components(G),
                desc="Reducing graph ",
                disable=not verbose,
                bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
            ):
                nodes = sorted(list(subg))
                # just assign the value in a single cell, because this will be representative of the clone
                tmp = nodes[0]
                for n in nodes:
                    nodes_names[n] = tmp  # keep so i can reference later
                if len(nodes) > 1:
                    G_ = G.subgraph(nodes).copy()
                    remove_edges[tmp] = [
                        (e[0], e[1])
                        for e in G_.edges(data=True)
                        if e[2]["weight"] > 0
                    ]
                    if len(remove_edges[tmp]) > 0:
                        G_.remove_edges_from(remove_edges[tmp])
                        for connected in nx.connected_components(G_):
                            vertexsizes[tmp].append(len(connected))
                        vertexsizes[tmp] = sorted(
                            vertexsizes[tmp], reverse=True
                        )
                    else:
                        vertexsizes[tmp] = [
                            1 for i in range(len(G_.edges(data=True)))
                        ]
                    if network_clustersize:
                        clustersizes[tmp] = len(vertexsizes[tmp])
                    else:
                        clustersizes[tmp] = len(nodes)
                else:
                    vertexsizes[tmp] = [1]
                    clustersizes[tmp] = [1]

            return (nodes_names, vertexsizes, clustersizes)
    else:
        raise TypeError("Input object must be of {}".format(Dandelion))


def diversity_gini(
    vdj_data: Dandelion | AnnData,
    groupby: str,
    metric: str | None = None,
    clone_key: str | None = None,
    min_size: int | None = None,
    n_boot: int = 1000,
    expanded_only: bool = False,
    use_contracted: bool = False,
    verbose: bool = False,
    **kwargs,
) -> tuple[pd.DataFrame, dict[list[float]], dict[list[float]]]:
    """
    Compute clones Gini indices.

    Parameters
    ----------
    vdj_data : Dandelion | AnnData
        Dandelion or AnnData object.
    groupby : str
        Column name to calculate the Gini indices on, for e.g. sample id, patient etc.
    metric : str | None, optional
        Metric to use for calculating Gini indices of clones.
        Accepts one of ['clone_network', 'clone_degree', 'clone_centrality'].
        Defaults to 'clone_centrality'.
    clone_key : str | None, optional
        Column name specifying the clone_id column in metadata.
    min_size : int | None, optional
        Minimum cell numbers to keep for diversity calculation. If None, defaults to size of smallest sample.
        Beware that this may lead to very small sample sizes and unreliable estimates if left as None.
    n_boot : int, optional
        Bootstrap iterations for calculations. Default is 1000.
    expanded_only : bool, optional
        Whether or not to calculate gini indices using expanded clones only. Default is False i.e. use all cells/clones.
    use_contracted : bool, optional
        Whether or not to perform the gini calculation after contraction of clone network.
        Only applies to calculation of clone size gini index. Default is False.
        This is to try and preserve the single-cell properties of the network.

    verbose : bool, optional
        whether to print progress.
    **kwargs
        passed to dandelion.tl.generate_network

    Returns
    -------
    tuple[pd.DataFrame, dict[list[float]]]
        pandas DataFrame holding summarised diversity estimation and the raw bootstrap results.

    Raises
    ------
    TypeError
        if not Dandelion class.
    ValueError
        if columns names don't exist.
    """
    logg.info("Calculating Gini indices")

    res, cluster_raw, vertex_raw = gini_indices(
        vdj_data,
        groupby=groupby,
        clone_key=clone_key,
        metric=metric,
        min_size=min_size,
        n_boot=n_boot,
        expanded_only=expanded_only,
        contracted=use_contracted,
        verbose=verbose,
        **kwargs,
    )

    return res, {"cluster_raw": cluster_raw, "vertex_raw": vertex_raw}


def chooseln(N, k) -> float:
    """
    R's lchoose in python
    from https://stackoverflow.com/questions/21767690/python-log-n-choose-k
    """
    return gammaln(N + 1) - gammaln(N - k + 1) - gammaln(k + 1)


def rarefun(y, sample_size) -> float:
    """
    Adapted from rarefun from vegan:
    https://github.com/vegandevs/vegan/blob/master/R/rarefy.R
    """
    res = []
    y = y[y > 0]
    J = np.sum(y)
    ldiv = chooseln(J, sample_size)
    for d in J - y:
        if d < sample_size:
            res.append(0)
        else:
            res.append(np.exp(chooseln(d, sample_size) - ldiv))
    out = np.sum(1 - np.array(res))
    return out


def drop_nan_values(df: pd.DataFrame) -> None:
    """
    Drop NaN values from a pandas DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame from which to drop NaN values.
    """
    if "nan" in df.index or np.nan in df.index:
        try:
            df.drop("nan", inplace=True)
        except:
            df.drop(np.nan, inplace=True)


def safe_bootstrap_summary(
    data_dict: dict,
    key: str,
    values: list[float],
    ci: float = 95,
) -> None:
    """
    Safely compute mean, standard deviation, and confidence interval
    from bootstrap values and update the dictionary.

    Parameters
    ----------
    data_dict : dict
        Dictionary to update.
    key : str
        Key identifying the group or sample.
    values : list of float
        Bootstrap values.
    ci : float, default=95
        Confidence interval width (percentile-based).
    """
    values = np.array(values, dtype=float)
    values = values[np.isfinite(values)]  # remove NaN/inf

    if len(values) == 0:
        mean, std, lower, upper = (np.nan, np.nan, np.nan, np.nan)
    else:
        mean = np.mean(values)
        std = np.std(values, ddof=1)
        alpha = (100 - ci) / 2
        lower, upper = np.percentile(values, [alpha, 100 - alpha])

    data_dict[key] = {
        "mean": mean,
        "std": std,
        f"lower_{ci}": lower,
        f"upper_{ci}": upper,
    }


def process_clone_network_stats(
    ddl_dat: Dandelion, expanded_only: bool, contracted: bool, verbose: bool
) -> tuple[dict, dict, dict]:
    """
    Process clone network statistics and calculate Gini indices.

    Parameters
    ----------
    ddl_dat : Any
        The Dandelion object to process.
    expanded_only : bool
        Whether to consider only expanded clones.
    contracted : bool
        Whether to use contracted network clusters.
    verbose : bool
        Whether to display progress information.

    Returns
    -------
    tuple[dict, dict, dict]
        Tuple containing dictionaries for node names, vertex sizes, and cluster sizes.
    """
    n_n, v_s, c_s = clone_networkstats(
        ddl_dat,
        expanded_only=expanded_only,
        network_clustersize=contracted,
        verbose=verbose,
    )
    g_c_v = defaultdict(dict)
    g_c_v_res, g_c_c_res = {}, {}

    for vs in v_s:
        v_sizes = np.array(v_s[vs])
        if len(v_sizes) > 1:
            v_sizes = np.append(v_sizes, 0)
        g_c_v[vs] = calculate_gini_index(v_sizes)
        for cell in n_n:
            g_c_v_res.update({cell: g_c_v[n_n[cell]]})

    c_sizes = np.array(sorted(list(flatten(c_s.values())), reverse=True))
    if len(c_sizes) > 1:
        c_sizes = np.append(c_sizes, 0)
    g_c_c = calculate_gini_index(c_sizes)
    for cell in n_n:
        g_c_c_res.update({cell: g_c_c})

    return g_c_v_res, g_c_c_res


def gini_indices(
    data: Dandelion,
    groupby: str,
    metric: str | None = None,
    clone_key: str | None = None,
    min_size: int | None = None,
    n_boot: int = 1000,
    expanded_only: bool = False,
    contracted: bool = False,
    verbose: bool = False,
    **kwargs,
) -> tuple[pd.DataFrame, dict[list[float]], dict[list[float]]]:
    """Gini indices."""
    if isinstance(data, AnnData):
        raise TypeError("Only Dandelion class object accepted.")
    clonekey = clone_key if clone_key is not None else "clone_id"
    met = metric if metric is not None else "clone_network"

    # --- Prepare data and filter groups
    data, _, groups, min_size = _prepare_diversity_groups(
        data, groupby, min_size
    )
    # filter the vdj object as well
    if isinstance(data, Dandelion):
        data = data[data.metadata[groupby].isin(groups)].copy()

    res1, res2, cluster_raw, vertex_raw = {}, {}, {}, {}
    for g in groups:
        # clone size distribution
        ddl_dat = data[data.metadata[groupby] == g].copy()

        sizelist, graphlist = [], []
        for _ in tqdm(
            range(0, n_boot),
            bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
            disable=not verbose,
        ):
            resample_sized = generate_network(
                ddl_dat,
                clone_key=clonekey,
                sample=min_size,
                force_replace=True,
                verbose=verbose,
                compute_layout=False,
                **kwargs,
            )
            if met == "clone_network":
                g_c_v_res, g_c_c_res = process_clone_network_stats(
                    resample_sized,
                    expanded_only=expanded_only,
                    contracted=contracted,
                    verbose=verbose,
                )
                resample_sized.metadata["clone_network_vertex_size_gini"] = (
                    pd.Series(g_c_v_res)
                )
                resample_sized.metadata["clone_network_cluster_size_gini"] = (
                    pd.Series(g_c_c_res)
                )
            elif met == "clone_centrality":
                clone_centrality(resample_sized)
            elif met == "clone_degree":
                clone_degree(resample_sized)
            else:
                raise ValueError(
                    "Unknown metric for calculating network stats. Please specify "
                    + "one of `clone_network`, `clone_centrality` or `clone_degree`."
                )
            # clone size gini
            _dat = resample_sized.metadata.copy()
            _tab = _dat[clonekey].value_counts()
            drop_nan_values(_tab)
            if met == "clone_network":
                sizelist.append(_dat[met + "_cluster_size_gini"].mean())
            else:
                clonesizecounts = np.array(_tab)
                clonesizecounts = clonesizecounts[clonesizecounts > 0]
                if len(clonesizecounts) > 1:
                    # append a single zero for lorenz curve calculation
                    clonesizecounts = np.append(clonesizecounts, 0)
                g_c = (
                    calculate_gini_index(clonesizecounts)
                    if len(clonesizecounts) > 0
                    else 0
                )
                sizelist.append(g_c)
            if met == "clone_network":
                graphlist.append(_dat[met + "_vertex_size_gini"].mean())
            else:
                # vertex closeness centrality or weighted degree distribution
                # only calculate for expanded clones. If including non-expanded clones, the centrality is just zero which doesn't help.
                connectednodes = resample_sized.metadata[met][
                    resample_sized.metadata[met] > 0
                ]
                graphcounts = np.array(connectednodes.value_counts())
                # graphcounts = np.append(graphcounts, 0) # if I add a  zero here, it will skew the results when the centrality measure is uniform.... so leave it out for now.
                g_c = (
                    calculate_gini_index(graphcounts)
                    if len(graphcounts) > 0
                    else 0
                )
                graphlist.append(g_c)
        cluster_raw[g] = sizelist.copy()
        vertex_raw[g] = graphlist.copy()
        safe_bootstrap_summary(res1, g, sizelist)
        safe_bootstrap_summary(res2, g, graphlist)
    res_df1 = pd.DataFrame.from_dict(res1).T
    res_df2 = pd.DataFrame.from_dict(res2).T
    res_df1.columns = [f"{met}_cluster_size_gini_" + x for x in res_df1.columns]
    res_df2.columns = [f"{met}_vertex_size_gini_" + x for x in res_df2.columns]
    res_df = pd.concat([res_df1, res_df2], axis=1)

    return res_df, cluster_raw, vertex_raw


def calculate_gini_index(
    values: np.ndarray, method: str = "trapezoids"
) -> float:
    """
    Calculate the Gini index for a given array of values.

    Parameters
    ----------
    values : np.ndarray
        Array of values to calculate the Gini index for.
    method : str, optional
        Method to use for Gini index calculation, by default "trapezoids".

    Returns
    -------
    float
        The calculated Gini index, or 0 if the index is negative or NaN.
    """
    gini = gini_index(values, method=method)
    return 0 if gini < 0 or np.isnan(gini) else gini


def calculate_chao1(values: np.ndarray) -> float:
    """
    Calculate the Chao1 estimate for a given array of values

    Parameters
    ----------
    values : np.ndarray
        Array of values to calculate the Chao1 estimates for.

    Returns
    -------
    float
        The calculated Chao1 estimates, or 0 if the index is negative or NaN.
    """
    chao1e = chao1(values)
    return 0 if chao1e < 0 or np.isnan(chao1e) else chao1e


def calculate_shannon_entropy(values: np.ndarray, normalize: bool) -> float:
    """
    Calculate the Shannon entropy for a given array of values

    Parameters
    ----------
    values : np.ndarray
        Array of values to calculate the Shannon entropys for.
    normalize : bool, optional
        Whether or not to return normalized Shannon Entropy according to https://math.stackexchange.com/a/945172. Default is True.
    Returns
    -------
    float
        The calculated Shannon entropys, or 0 if the index is negative or NaN.
    """
    if normalize:
        if len(values) == 1:
            return 0
        else:
            values_freqs = values / np.sum(values)
            return -np.sum(
                (values_freqs * np.log(values_freqs))
                / np.log(len(values_freqs))
            )
    else:
        return shannon(values)


def rename_result_column(
    res_df: pd.DataFrame,
    diversity_mode: str,
) -> None:
    """Processes the output result"""
    res_df.columns = [
        f"clone_size_{diversity_mode}_{x}" for x in res_df.columns
    ]


def min_sample_size(
    df: pd.DataFrame,
    col: str,
):
    """
    Calculate the minimum group size for downsampling and log related information.

    Parameters
    ----------
    _df : pd.DataFrame
        Metadata DataFrame containing the grouping information.
    col : str
        Column name in `_df` to group by.

    Returns
    -------
    int
        The determined minimum group size (`min_size`).
    """
    # Determine minimum group size
    min_size = df[col].value_counts().min()
    return min_size


def diversity_estimates(
    vdj_data: Dandelion | AnnData,
    groupby: str,
    method: Literal["chao1", "shannon", "gini"] = "chao1",
    clone_key: str | None = None,
    normalize: bool = True,
    min_size: int | None = None,
    n_boot: int = 1000,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Compute clones Chao1 estimates.

    Parameters
    ----------
    vdj_data : Dandelion | AnnData
        Dandelionr AnnData object.
    groupby : str
        Column name to calculate the Chao1 estimates on, for e.g. sample id, patient etc.
    method : Literal["chao1", "shannon", "gini"], optional
        Diversity metric to compute.
    clone_key : str | None, optional
        Column name specifying the clone_id column in metadata
    normalize : bool, optional
        Whether or not to return normalized Shannon Entropy according to https://math.stackexchange.com/a/945172. Default is True.
    min_size : int | None, optional
        Minimum cell numbers to keep for diversity calculation. If None, defaults to size of smallest sample.
        Beware that this may lead to very small sample sizes and unreliable estimates if left as None.
    n_boot : int, optional
        Bootstrap iterations for calculations. Default is 1000.

    verbose : bool, optional
        whether to print progress.

    Returns
    -------
    pd.DataFrame
        pandas DataFrame holding diversity information.
    """
    res, res_raw = estimate_diversity(
        vdj_data,
        groupby=groupby,
        clone_key=clone_key,
        metric=method,
        normalize=normalize,
        min_size=min_size,
        n_boot=n_boot,
        verbose=verbose,
    )

    return res, res_raw


def estimate_diversity(
    data: Dandelion | AnnData,
    groupby: str,
    clone_key: str | None = None,
    metric: Literal["chao1", "shannon", "gini"] = "chao1",
    normalize: bool = True,  # used only if metric == "shannon"
    min_size: int | None = None,
    n_boot: int = 1000,
    verbose: bool = False,
) -> tuple[pd.DataFrame, dict[list[float]]]:
    """
    Estimate immune repertoire diversity metrics (Chao1 or Shannon).

    Parameters
    ----------
    data : Dandelion | AnnData
        Input repertoire or annotated single-cell data.
    groupby : str
        Column name used to group cells/samples for diversity estimation.
    clone_key : str, optional
        Column containing clone identifiers. Defaults to "clone_id".
    metric : {"chao1", "shannon"}, default="chao1"
        Diversity metric to compute.
    normalize : bool, default=True
        If True, returns normalized Shannon entropy (ignored for Chao1).
    min_size : int, optional
        Minimum group size required for diversity calculation.
    n_boot : int, default=1000
        Number of bootstrap iterations.
    verbose : bool, default=False
        If True, show progress bars.

    Returns
    -------
    tuple[pd.DataFrame, dict[list[float]]]
        pandas DataFrame holding summarised diversity estimation and the raw bootstrap results.

    """
    clonekey = clone_key if clone_key is not None else "clone_id"

    # Determine diversity mode label
    if metric.lower() == "shannon":
        diversity_mode = "normalized_shannon" if normalize else "shannon"
    else:
        diversity_mode = metric.lower()

    # --- Prepare data and filter groups
    data, _metadata, groups, min_size = _prepare_diversity_groups(
        data, groupby, min_size
    )

    # --- Subset the actual data
    if isinstance(data, Dandelion):
        data = data[data.metadata[groupby].isin(groups)].copy()
    elif isinstance(data, AnnData):
        data = data[data.obs[groupby].isin(groups)].copy()

    # --- Compute diversity estimates via bootstrapping
    res, res_raw = {}, {}
    for g in groups:
        _dat = _metadata[_metadata[groupby] == g]
        values = []

        for _ in tqdm(
            range(n_boot),
            bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
            disable=not verbose,
        ):
            _sample = _dat.sample(min_size, replace=True)
            _tab = _sample[clonekey].value_counts()
            drop_nan_values(_tab)
            clone_sizes = np.array(_tab)
            clone_sizes = clone_sizes[clone_sizes > 0]

            if len(clone_sizes) == 0:
                val = 0
            elif metric.lower() == "chao1":
                val = calculate_chao1(clone_sizes)
            elif metric.lower() == "shannon":
                val = calculate_shannon_entropy(clone_sizes, normalize)
            elif metric.lower() == "gini":
                val = calculate_gini_index(clone_sizes)
            values.append(val)

        res_raw[g] = values.copy()
        safe_bootstrap_summary(res, g, values)

    # --- Format result
    res_df = pd.DataFrame.from_dict([res])
    rename_result_column(res_df, diversity_mode)

    return res_df, res_raw


def _prepare_diversity_groups(
    data: Dandelion | AnnData,
    groupby: str,
    min_size: int | None = None,
):
    """Shared logic for preparing metadata and filtering valid groups."""
    if isinstance(data, AnnData):
        _metadata = data.obs.copy()
    elif isinstance(data, Dandelion):
        _metadata = data.metadata.copy()
    else:
        raise TypeError("data must be an AnnData or Dandelion object")

    _metadata[groupby] = _metadata[groupby].astype("category")
    _metadata[groupby] = _metadata[groupby].cat.remove_unused_categories()
    groups = list(set(_metadata[groupby]))

    if min_size is None:
        min_size = min_sample_size(df=_metadata, col=groupby)
    group_counts = _metadata[groupby].value_counts()
    valid_groups = group_counts[group_counts >= min_size].index.tolist()
    if len(valid_groups) < 1:
        raise ValueError(
            f"No groups have sufficient size for diversity calculation. "
            f"Please choose a smaller min_size than {min_size}."
        )
    if len(valid_groups) < len(groups):
        logg.warning(
            f"The following groups were excluded (<{min_size} cells): "
            + ", ".join(set(groups) - set(valid_groups))
        )
    groups = valid_groups
    _metadata = _metadata[_metadata[groupby].isin(groups)].copy()

    return data, _metadata, groups, min_size
