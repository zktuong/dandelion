#!/usr/bin/env python
import numpy as np
import networkx as nx
import pandas as pd
import warnings

from anndata import AnnData
from time import sleep
from tqdm import tqdm
from scanpy import logging as logg
from scipy.special import gammaln
from typing import Literal

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
    diversity_key: str | None = None,
    verbose: bool = False,
) -> dict:
    """
    Return rarefaction predictions for cell numbers vs clone size.

    Parameters
    ----------
    vdj_data : Dandelion | AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to split the calculation of clone numbers for a given number of cells for e.g. sample, patient etc.
    clone_key : str | None, optional
        Column name specifying the clone_id column in metadata/obs.
    diversity_key : str | None, optional
        key for 'diversity' results in AnnData's `.uns`.
    verbose : bool, optional
        whether to print progress.

    Returns
    -------
    dict
        Rarefaction predictions for cell numbers vs clone size.
    """
    start = logg.info("Constructing rarefaction curve")

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

    diversitykey = diversity_key if diversity_key is not None else "diversity"

    if isinstance(vdj_data, AnnData):
        if diversitykey not in vdj_data.uns:
            vdj_data.uns[diversitykey] = {}
        vdj_data.uns[diversitykey] = {
            "rarefaction_cells_x": pred,
            "rarefaction_clones_y": y,
        }
    logg.info(
        " finished",
        time=start,
        deep=("updated `.uns` with rarefaction curves.\n"),
    )
    if isinstance(vdj_data, Dandelion):
        return {"rarefaction_cells_x": pred, "rarefaction_clones_y": y}


def clone_diversity(
    vdj_data: Dandelion | AnnData,
    groupby: str,
    method: Literal["gini", "chao1", "shannon"] = "gini",
    metric: Literal["clone_network", "clone_degree", "clone_centrality"] = None,
    clone_key: str | None = None,
    return_table: bool = False,
    diversity_key: str | None = None,
    resample: bool = False,
    downsample: int | None = None,
    n_resample: int = 50,
    normalize: bool = True,
    reconstruct_network: bool = True,
    expanded_only: bool = False,
    use_contracted: bool = False,
    key_added: str | None = None,
    verbose: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """
    Compute B cell clones diversity : Gini indices, Chao1 estimates, or Shannon entropy.

    Parameters
    ----------
    vdj_data : Dandelion | AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the gini indices on, for e.g. sample, patient etc.
    method : Literal["gini", "chao1", "shannon"], optional
        Method for diversity estimation. Either one of ['gini', 'chao1', 'shannon'].
    metric : Literal["clone_network", "clone_degree", "clone_centrality"], optional
        Metric to use for calculating Gini indices of clones.
        Accepts one of ['clone_network', 'clone_degree', 'clone_centrality'].
        `None` defaults to 'clone_network'.
    clone_key : str | None, optional
        Column name specifying the clone_id column in metadata.
    return_table : bool, optional
        If True, a `pandas` data frame is returned.
        If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str | None, optional
        key for 'diversity' results in `.uns`.
    resample : bool, optional
        Whether or not to randomly sample cells without replacement to
        the minimum size of groups for the diversity calculation. Default is False.
    downsample : int | None, optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    n_resample : int, optional
        Number of times to perform resampling. Default is 50.
    normalize : bool, optional
        Whether or not to return normalized Shannon Entropy according to https://math.stackexchange.com/a/945172.
        Default is True.
    reconstruct_network : bool, optional
        Whether or not to reconstruct the network for Gini Index based measures.
        Default is True and will reconstruct for each group specified by groupby option.
    expanded_only : bool, optional
        Whether or not to calculate gini indices using expanded clones only. Default is False i.e. use all cells/clones.
    use_contracted : bool, optional
        Whether or not to perform the gini calculation after contraction of clone network.
        Only applies to calculation of clone size gini index. Default is False.
        This is to try and preserve the single-cell properties of the network.
    key_added : str | None, optional
        column names for output.
    verbose : bool, optional
        whether to print progress.
    **kwargs
        passed to dandelion.tl.generate_network

    Returns
    -------
    pd.DataFrame
        Pandas DataFrame holding diversity information.
    """
    if downsample is not None:
        resample = True

    if method == "gini":
        if return_table:
            return diversity_gini(
                vdj_data,
                groupby=groupby,
                metric=metric,
                clone_key=clone_key,
                return_table=return_table,
                diversity_key=diversity_key,
                resample=resample,
                n_resample=n_resample,
                downsample=downsample,
                reconstruct_network=reconstruct_network,
                expanded_only=expanded_only,
                use_contracted=use_contracted,
                key_added=key_added,
                verbose=verbose,
                **kwargs,
            )
        else:
            diversity_gini(
                vdj_data,
                groupby=groupby,
                metric=metric,
                clone_key=clone_key,
                return_table=return_table,
                diversity_key=diversity_key,
                resample=resample,
                n_resample=n_resample,
                downsample=downsample,
                reconstruct_network=reconstruct_network,
                expanded_only=expanded_only,
                use_contracted=use_contracted,
                key_added=key_added,
                verbose=verbose,
                **kwargs,
            )
    if method == "chao1":
        if return_table:
            return diversity_chao1(
                vdj_data,
                groupby=groupby,
                clone_key=clone_key,
                return_table=return_table,
                diversity_key=diversity_key,
                resample=resample,
                n_resample=n_resample,
                downsample=downsample,
                key_added=key_added,
                verbose=verbose,
            )
        else:
            diversity_chao1(
                vdj_data,
                groupby=groupby,
                clone_key=clone_key,
                return_table=return_table,
                diversity_key=diversity_key,
                resample=resample,
                n_resample=n_resample,
                downsample=downsample,
                key_added=key_added,
                verbose=verbose,
            )
    if method == "shannon":
        if return_table:
            return diversity_shannon(
                vdj_data,
                groupby=groupby,
                clone_key=clone_key,
                return_table=return_table,
                diversity_key=diversity_key,
                resample=resample,
                n_resample=n_resample,
                normalize=normalize,
                downsample=downsample,
                key_added=key_added,
                verbose=verbose,
            )
        else:
            diversity_shannon(
                vdj_data,
                groupby=groupby,
                clone_key=clone_key,
                return_table=return_table,
                diversity_key=diversity_key,
                resample=resample,
                n_resample=n_resample,
                normalize=normalize,
                downsample=downsample,
                key_added=key_added,
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
    start = logg.info("Calculating vertex size of nodes after contraction")

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
    return_table: bool = False,
    resample: bool = False,
    n_resample: int = 50,
    downsample: int | None = None,
    reconstruct_network: bool = True,
    expanded_only: bool = False,
    use_contracted: bool = False,
    key_added: str | None = None,
    verbose: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """
    Compute clones Gini indices.

    Parameters
    ----------
    vdj_data : Dandelion | AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Gini indices on, for e.g. sample, patient etc.
    metric : str | None, optional
        Metric to use for calculating Gini indices of clones.
        Accepts one of ['clone_network', 'clone_degree', 'clone_centrality'].
        Defaults to 'clone_centrality'.
    clone_key : str | None, optional
        Column name specifying the clone_id column in metadata.
    return_table : bool, optional
        If True, a `pandas` data frame is returned.
        If False, function will try to populate the input object's metadata/obs slot.
    resample : bool, optional
        Whether or not to randomly sample cells without replacement to
        the minimum size of groups for the diversity calculation.
        Default is False. Resampling will automatically trigger reconstruction of network.
    n_resample : int, optional
        Number of times to perform resampling. Default is 50.
    downsample : int | None, optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    reconstruct_network : bool, optional
        Whether or not to reconstruct the network for Gini Index based measures.
        Default is True and will reconstruct for each group specified by groupby option.
    expanded_only : bool, optional
        Whether or not to calculate gini indices using expanded clones only. Default is False i.e. use all cells/clones.
    use_contracted : bool, optional
        Whether or not to perform the gini calculation after contraction of clone network.
        Only applies to calculation of clone size gini index. Default is False.
        This is to try and preserve the single-cell properties of the network.
    key_added : str | None, optional
        column names for output.
    verbose : bool, optional
        whether to print progress.
    **kwargs
        passed to dandelion.tl.generate_network

    Returns
    -------
    pd.DataFrame
        pandas DataFrame holding diversity information.

    Raises
    ------
    TypeError
        if not Dandelion class.
    ValueError
        if columns names don't exist.
    """
    start = logg.info("Calculating Gini indices")

    res = gini_indices(
        vdj_data,
        groupby=groupby,
        clone_key=clone_key,
        metric=metric,
        resample=resample,
        n_resample=n_resample,
        downsample=downsample,
        reconstruct_network=reconstruct_network,
        expanded_only=expanded_only,
        contracted=use_contracted,
        key_added=key_added,
        verbose=verbose,
        **kwargs,
    )

    if return_table:
        logg.info(" finished", time=start)
        return res
    else:
        transfer_diversity_results(vdj_data, res, groupby)
        if isinstance(vdj_data, Dandelion):
            logg.info(
                " finished",
                time=start,
                deep=("updated `.metadata` with Gini indices.\n"),
            )


def diversity_chao1(
    vdj_data: Dandelion | AnnData,
    groupby: str,
    clone_key: str | None = None,
    return_table: bool = False,
    diversity_key: str | None = None,
    resample: bool = False,
    n_resample: int = 50,
    downsample: int | None = None,
    key_added: str | None = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Compute clones Chao1 estimates.

    Parameters
    ----------
    vdj_data : Dandelion | AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Chao1 estimates on, for e.g. sample, patient etc.
    clone_key : str | None, optional
        Column name specifying the clone_id column in metadata.
    return_table : bool, optional
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's
        metadata/obs slot.
    diversity_key : str | None, optional
        key for 'diversity' results in `.uns`.
    resample : bool, optional
        Whether or not to randomly sample cells without replacement to the minimum size of groups for the diversity
        calculation. Default is False.
    n_resample : int, optional
        Number of times to perform resampling. Default is 50.
    downsample : int | None, optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    key_added : str | None, optional
        column names for output.
    verbose : bool, optional
        whether to print progress.

    Returns
    -------
    pd.DataFrame
        pandas DataFrame holding diversity information.
    """
    start = logg.info("Calculating Chao1 estimates")

    res = chao1_estimates(
        vdj_data,
        groupby=groupby,
        clone_key=clone_key,
        resample=resample,
        n_resample=n_resample,
        downsample=downsample,
        key_added=key_added,
        verbose=verbose,
    )

    diversitykey = diversity_key if diversity_key is not None else "diversity"

    if isinstance(vdj_data, AnnData):
        if diversitykey not in vdj_data.uns:
            vdj_data.uns[diversitykey] = {}
        vdj_data.uns[diversitykey].update({"chao1": res})

    if return_table:
        if isinstance(vdj_data, AnnData):
            logg.info(
                " finished",
                time=start,
                deep=("updated `.uns` with Chao1 estimates.\n"),
            )
        else:
            logg.info(" finished", time=start)
        return res
    else:
        transfer_diversity_results(vdj_data, res, groupby)
        if isinstance(vdj_data, Dandelion):
            logg.info(
                " finished",
                time=start,
                deep=("updated `.metadata` with Chao1 estimates.\n"),
            )
        elif isinstance(vdj_data, AnnData):
            logg.info(
                " finished",
                time=start,
                deep=("updated `.obs` and `.uns` with Chao1 estimates.\n"),
            )


def diversity_shannon(
    vdj_data: Dandelion | AnnData,
    groupby: str,
    clone_key: str | None = None,
    return_table: bool = False,
    diversity_key: str | None = None,
    resample: bool = False,
    n_resample: int = 50,
    normalize: bool = True,
    key_added: str | None = None,
    downsample: int | None = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Compute clones Shannon entropy.

    Parameters
    ----------
    vdj_data : Dandelion | AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Shannon entropy on, for e.g. sample, patient etc.
    clone_key : str | None, optional
        Column name specifying the clone_id column in metadata.
    return_table : bool, optional
        If True, a `pandas` data frame is returned. If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str | None, optional
        key for 'diversity' results in `.uns`.
    resample : bool, optional
        Whether or not to randomly sample cells without replacement to the minimum size of groups for the diversity calculation. Default is False.
    n_resample : int, optional
        Number of times to perform resampling. Default is 50.
    normalize : bool, optional
        Whether or not to return normalized Shannon Entropy according to https://math.stackexchange.com/a/945172. Default is True.
    key_added : str | None, optional
        column names for output.
    downsample : int | None, optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    verbose : bool, optional
        whether to print progress.

    Returns
    -------
    pd.DataFrame
        pandas DataFrame holding diversity information.
    """
    start = logg.info("Calculating Shannon entropy")

    res = shannon_entropy(
        vdj_data,
        groupby=groupby,
        clone_key=clone_key,
        resample=resample,
        n_resample=n_resample,
        normalize=normalize,
        key_added=key_added,
        downsample=downsample,
        verbose=verbose,
    )

    diversitykey = diversity_key if diversity_key is not None else "diversity"

    if isinstance(vdj_data, AnnData):
        if diversitykey not in vdj_data.uns:
            vdj_data.uns[diversitykey] = {}
        vdj_data.uns[diversitykey].update({"shannon": res})

    if return_table:
        if isinstance(vdj_data, AnnData):
            if normalize:
                logg.info(
                    " finished",
                    time=start,
                    deep=("updated `.uns` with normalized Shannon entropy.\n"),
                )
            else:
                logg.info(
                    " finished",
                    time=start,
                    deep=("updated `.uns` with Shannon entropy.\n"),
                )
        else:
            logg.info(" finished", time=start)
        return res
    else:
        transfer_diversity_results(vdj_data, res, groupby)
        if isinstance(vdj_data, Dandelion):
            if normalize:
                logg.info(
                    " finished",
                    time=start,
                    deep=(
                        "updated `.metadata` with normalized Shannon entropy.\n"
                    ),
                )
            else:
                logg.info(
                    " finished",
                    time=start,
                    deep=("updated `.metadata` with Shannon entropy.\n"),
                )
        elif isinstance(vdj_data, AnnData):
            if normalize:
                logg.info(
                    " finished",
                    time=start,
                    deep=(
                        "updated `.obs` and `.uns` with normalized Shannon entropy.\n"
                    ),
                )
            else:
                logg.info(
                    " finished",
                    time=start,
                    deep=("updated `.obs` and `.uns` with Shannon entropy.\n"),
                )


def chooseln(N, k) -> float:
    """
    R's lchoose in python
    from https://stackoverflow.com/questions/21767690/python-log-n-choose-k
    """
    return gammaln(N + 1) - gammaln(N - k + 1) - gammaln(k + 1)


def rarefun(y, sample) -> float:
    """
    Adapted from rarefun from vegan:
    https://github.com/vegandevs/vegan/blob/master/R/rarefy.R
    """
    res = []
    y = y[y > 0]
    J = np.sum(y)
    ldiv = chooseln(J, sample)
    for d in J - y:
        if d < sample:
            res.append(0)
        else:
            res.append(np.exp(chooseln(d, sample) - ldiv))
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


def safe_average_update(data_dict: dict, key: str, values: list[float]) -> None:
    """
    Safely calculate the average of a list of values and update the dictionary.

    Parameters
    ----------
    data_dict : dict
        The dictionary to update.
    key : Any
        The key to update in the dictionary.
    values : list of float
        The list of values to average.
    """
    try:
        average = sum(values) / len(values)
    except ZeroDivisionError:
        average = 0
    data_dict.update({key: average})


def process_clone_network_stats(
    ddl_dat: Dandelion, expanded_only: bool, contracted: bool, verbose: bool
) -> tuple[dict, dict, dict]:
    """
    Process clone network statistics and calculate Gini indices.

    Parameters
    ----------
    ddl_dat : Any
        The `Dandelion` object to process.
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
    resample: bool = False,
    n_resample: int = 50,
    downsample: int | None = None,
    reconstruct_network: bool = True,
    expanded_only: bool = False,
    contracted: bool = False,
    key_added: str | None = None,
    verbose: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """Gini indices."""
    if isinstance(data, AnnData):
        raise TypeError("Only Dandelion class object accepted.")
    elif isinstance(data, Dandelion):
        _metadata = data.metadata.copy()
    clonekey = clone_key if clone_key is not None else "clone_id"
    met = metric if metric is not None else "clone_network"
    # split up the table by groupby
    _metadata[groupby] = _metadata[groupby].astype("category")
    _metadata[groupby] = _metadata[groupby].cat.remove_unused_categories()
    groups = list(set(_metadata[groupby]))

    minsize = calculate_minsize_and_log(
        df=_metadata, col=groupby, downsample=downsample, resample=resample
    )

    res1 = {}
    if met == "clone_network":
        logg.info(
            "Computing Gini indices for cluster and vertex size using network."
        )
        if not reconstruct_network:
            n_n, v_s, c_s = clone_networkstats(
                data,
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
            c_sizes = np.array(
                np.array(sorted(list(flatten(c_s.values())), reverse=True))
            )
            if len(c_sizes) > 1:
                c_sizes = np.append(c_sizes, 0)
            g_c_c = calculate_gini_index(c_sizes)
            for cell in n_n:
                g_c_c_res.update({cell: g_c_c})
            data.metadata["clone_network_vertex_size_gini"] = pd.Series(
                g_c_v_res
            )
            data.metadata["clone_network_cluster_size_gini"] = pd.Series(
                g_c_c_res
            )
    elif met == "clone_centrality":
        logg.info(
            "Computing gini indices for clone size using metadata and node closeness centrality using network."
        )
        clone_centrality(data)
    elif met == "clone_degree":
        logg.info(
            "Computing gini indices for clone size using metadata and node degree using network."
        )
        clone_degree(data)
    _metadata = data.metadata.copy()
    _data = data.data.copy()
    res2 = {}
    if resample:
        logg.info(
            "Downsampling each group specified in `{}` to {} cells for calculating gini indices.".format(
                groupby, minsize
            )
        )
    for g in groups:
        # clone size distribution
        _dat = _metadata[_metadata[groupby] == g]
        __data = _data[_data["cell_id"].isin(list(_dat.index))]
        ddl_dat = Dandelion(__data, metadata=_dat)
        if resample:
            sizelist = []
            if isinstance(data, Dandelion):
                graphlist = []
            for i in tqdm(
                range(0, n_resample),
                bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
                disable=not verbose,
            ):
                if isinstance(data, Dandelion):
                    resampled = generate_network(
                        ddl_dat,
                        clone_key=clonekey,
                        downsample=minsize,
                        verbose=verbose,
                        compute_layout=False,
                        **kwargs,
                    )
                    if met == "clone_network":
                        g_c_v_res, g_c_c_res = process_clone_network_stats(
                            resampled,
                            expanded_only=expanded_only,
                            contracted=contracted,
                            verbose=verbose,
                        )
                        resampled.metadata["clone_network_vertex_size_gini"] = (
                            pd.Series(g_c_v_res)
                        )
                        resampled.metadata[
                            "clone_network_cluster_size_gini"
                        ] = pd.Series(g_c_c_res)
                    elif met == "clone_centrality":
                        clone_centrality(resampled)
                    elif met == "clone_degree":
                        clone_degree(resampled)
                    else:
                        raise ValueError(
                            "Unknown metric for calculating network stats. Please specify "
                            + "one of `clone_network`, `clone_centrality` or `clone_degree`."
                        )
                    # clone size gini
                    _dat = resampled.metadata.copy()
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
                        connectednodes = resampled.metadata[met][
                            resampled.metadata[met] > 0
                        ]
                        graphcounts = np.array(connectednodes.value_counts())
                        # graphcounts = np.append(graphcounts, 0) # if I add a  zero here, it will skew the results when the centrality measure is uniform.... so leave it out for now.
                        g_c = (
                            calculate_gini_index(graphcounts)
                            if len(graphcounts) > 0
                            else 0
                        )
                        graphlist.append(g_c)

            safe_average_update(res1, g, sizelist)
            if "graphlist" in locals():
                safe_average_update(res2, g, graphlist)
        else:
            _tab = _dat[clonekey].value_counts()
            drop_nan_values(_tab)
            if met != "clone_network":
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
                res1.update({g: g_c})
            if isinstance(data, Dandelion):
                if met == "clone_network":
                    if reconstruct_network:
                        generate_network(
                            ddl_dat,
                            clone_key=clonekey,
                            verbose=verbose,
                            compute_layout=False,
                            **kwargs,
                        )
                        g_c_v_res, g_c_c_res = process_clone_network_stats(
                            ddl_dat,
                            expanded_only=expanded_only,
                            contracted=contracted,
                            verbose=verbose,
                        )
                        res2.update({g: pd.Series(g_c_v_res).mean()})
                        res1.update({g: pd.Series(g_c_c_res).mean()})
                    else:
                        res2.update({g: _dat[met + "_vertex_size_gini"].mean()})
                        res1.update(
                            {g: _dat[met + "_cluster_size_gini"].mean()}
                        )
                else:
                    # vertex closeness centrality or weighted degree distribution
                    # only calculate for expanded clones. If including non-expanded clones, the centrality is
                    # just zero which doesn't help.
                    connectednodes = _dat[met][_dat[met] > 0]
                    graphcounts = np.array(connectednodes.value_counts())
                    # graphcounts = np.append(graphcounts, 0) # if I add a  zero here, it will skew the results
                    # when the centrality measure is uniform.... so leave it out for now.
                    g_c = (
                        calculate_gini_index(graphcounts)
                        if len(graphcounts) > 0
                        else 0
                    )
                    res2.update({g: g_c})
    if "res2" in locals():
        res_df = pd.DataFrame.from_dict([res1, res2]).T
        if key_added is None:
            if met == "clone_network":
                res_df.columns = [
                    met + "_cluster_size_gini",
                    met + "_vertex_size_gini",
                ]
            else:
                res_df.columns = ["clone_size_gini", met + "_gini"]
        else:
            if not type(key_added) is list:
                key_added = [key_added]
            if len(key_added) == len(res_df.columns):
                res_df.columns = key_added
            else:
                raise ValueError(
                    "Please provide {} key(s) for new column names.".format(
                        len(res_df.columns)
                    )
                )
    else:
        res_df = pd.DataFrame.from_dict([res1]).T
        if key_added is None:
            res_df.columns = ["clone_size_gini"]
        else:
            if not type(key_added) is list:
                key_added = [key_added]
            if len(key_added) == len(res_df.columns):
                res_df.columns = key_added
            else:
                raise ValueError(
                    "Please provide {} key(s) for new column names.".format(
                        len(res_df.columns)
                    )
                )
    return res_df


def chao1_estimates(
    data: Dandelion | AnnData,
    groupby: str,
    clone_key: str | None = None,
    resample: bool = False,
    n_resample: int = 50,
    downsample: int | None = None,
    key_added: str | None = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Chao1 estimates."""
    if isinstance(data, AnnData):
        _metadata = data.obs.copy()
    elif isinstance(data, Dandelion):
        _metadata = data.metadata.copy()
    clonekey = clone_key if clone_key is not None else "clone_id"
    diversity_mode = "chao1"
    # split up the table by groupby
    _metadata[groupby] = _metadata[groupby].astype("category")
    _metadata[groupby] = _metadata[groupby].cat.remove_unused_categories()
    groups = list(set(_metadata[groupby]))

    minsize = calculate_minsize_and_log(
        df=_metadata, col=groupby, downsample=downsample, resample=resample
    )

    res1 = {}
    for g in groups:
        # clone size distribution
        _dat = _metadata[_metadata[groupby] == g]
        if resample:
            sizelist = []
            for i in tqdm(
                range(0, n_resample),
                bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
                disable=not verbose,
            ):
                _dat = _dat.sample(minsize)
                _tab = _dat[clonekey].value_counts()
                drop_nan_values(_tab)
                clonesizecounts = np.array(_tab)
                clonesizecounts = clonesizecounts[clonesizecounts > 0]
                g_c = (
                    calculate_chao1(clonesizecounts)
                    if len(clonesizecounts) > 0
                    else 0
                )
                sizelist.append(g_c)
            safe_average_update(res1, g, sizelist)
        else:
            _tab = _dat[clonekey].value_counts()
            drop_nan_values(_tab)
            clonesizecounts = np.array(_tab)
            clonesizecounts = clonesizecounts[clonesizecounts > 0]
            g_c = (
                calculate_chao1(clonesizecounts)
                if len(clonesizecounts) > 0
                else 0
            )
            res1.update({g: g_c})

    res_df = pd.DataFrame.from_dict([res1]).T
    rename_result_column(res_df, diversity_mode, key_added)

    return res_df


def shannon_entropy(
    data: Dandelion | AnnData,
    groupby: str,
    clone_key: str | None = None,
    resample: bool = False,
    n_resample: int = 50,
    normalize: bool = True,
    downsample: int | None = None,
    key_added: str | None = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Shannon entropy."""
    if isinstance(data, AnnData):
        _metadata = data.obs.copy()
    elif isinstance(data, Dandelion):
        _metadata = data.metadata.copy()
    clonekey = clone_key if clone_key is not None else "clone_id"
    diversity_mode = "shannon" if not normalize else "normalized_shannon"
    # split up the table by groupby
    _metadata[groupby] = _metadata[groupby].astype("category")
    _metadata[groupby] = _metadata[groupby].cat.remove_unused_categories()
    groups = list(set(_metadata[groupby]))

    minsize = calculate_minsize_and_log(
        df=_metadata, col=groupby, downsample=downsample, resample=resample
    )

    res1 = {}
    for g in groups:
        # clone size distribution
        _dat = _metadata[_metadata[groupby] == g]
        if resample:
            sizelist = []
            for i in tqdm(
                range(0, n_resample),
                bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
                disable=not verbose,
            ):
                _dat = _dat.sample(minsize)
                _tab = _dat[clonekey].value_counts()
                drop_nan_values(_tab)
                clonesizecounts = np.array(_tab)
                clonesizecounts = clonesizecounts[clonesizecounts > 0]
                g_c = (
                    calculate_shannon_entropy(clonesizecounts, normalize)
                    if len(clonesizecounts) > 0
                    else 0
                )
                sizelist.append(g_c)
            safe_average_update(res1, g, sizelist)
        else:
            _tab = _dat[clonekey].value_counts()
            drop_nan_values(_tab)
            clonesizecounts = np.array(_tab)
            clonesizecounts = clonesizecounts[clonesizecounts > 0]
            g_c = (
                calculate_shannon_entropy(clonesizecounts, normalize)
                if len(clonesizecounts) > 0
                else 0
            )
            res1.update({g: g_c})

    res_df = pd.DataFrame.from_dict([res1]).T
    rename_result_column(res_df, diversity_mode, key_added)

    return res_df


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


def transfer_diversity_results(
    vdj_data: Dandelion | AnnData, diversity_results: pd.DataFrame, groupby: str
) -> None:
    """Transfer diversity results."""
    if isinstance(vdj_data, AnnData):
        _metadata = vdj_data.obs.copy()
    elif isinstance(vdj_data, Dandelion):
        _metadata = vdj_data.metadata.copy()

    groups = list(set(_metadata[groupby]))
    for c in diversity_results.columns:
        _metadata[c] = np.nan
        for g in groups:
            for i in _metadata.index:
                if _metadata.at[i, groupby] == g:
                    _metadata.at[i, c] = diversity_results[c][g]
    if isinstance(vdj_data, AnnData):
        vdj_data.obs = _metadata.copy()
    elif isinstance(vdj_data, Dandelion):
        vdj_data.metadata = _metadata.copy()


def rename_result_column(
    res_df: pd.DataFrame,
    diversity_mode: str,
    key_added: list[str] | str | None = None,
) -> None:
    """Processes the output result"""
    if key_added is None:
        res_df.columns = ["clone_size_" + diversity_mode]
    else:
        if isinstance(key_added, list):
            res_df.columns = key_added[0]
        else:
            res_df.columns = [key_added]


def calculate_minsize_and_log(
    df: pd.DataFrame,
    col: str,
    downsample: int | None = None,
    resample: int | bool | None = False,
):
    """
    Calculate the minimum group size for downsampling and log related information.

    Parameters
    ----------
    _df : pd.DataFrame
        Metadata DataFrame containing the grouping information.
    col : str
        Column name in `_df` to group by.
    downsample : int | None, optional
        Downsampling size. If None or larger than the smallest group size, defaults to the smallest group size.
    resample : bool | int | None, optional
        Whether to log a message about resampling.
    logg : logging.Logger, optional
        Logger instance for logging messages. If None, no logs will be emitted.

    Returns
    -------
    int
        The determined minimum group size (`minsize`).

    Notes
    -----
    - If `downsample` is greater than the smallest group size, the function logs a warning and defaults
      to the smallest group size.
    - If `minsize` is less than 100, a cautionary message is logged about potential issues with
      diversity measures.
    - When `resample` is True, a message is logged about the downsampling process.
    """
    # Determine minimum group size
    if downsample is None:
        minsize = df[col].value_counts().min()
    else:
        minsize = downsample
        if minsize > df[col].value_counts().min():

            logg.info(
                "Downsampling size provided of {} was larger than the smallest group size. ".format(
                    downsample
                )
                + "Defaulting to the smallest group size for downsampling."
            )
            minsize = df[col].value_counts().min()

    # Log a warning if the minimum size is too small
    if minsize < 100:

        logg.info(
            "The minimum cell numbers when grouped by {} is {}.".format(
                col, minsize
            )
            + " Exercise caution when interpreting diversity measures."
        )

    # Log information about resampling
    if resample:

        logg.info(
            "Downsampling each group specified in `{}` to {} cells for calculating Shannon entropy.".format(
                col, minsize
            )
        )

    return minsize
