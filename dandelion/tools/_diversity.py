#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-08-13 21:08:53
# @Last Modified by:   Kelvin
# @Last Modified time: 2022-10-27 10:19:36
"""diversity module."""
import numpy as np
import networkx as nx
import pandas as pd
import warnings

from anndata import AnnData
from time import sleep
from tqdm import tqdm
from scanpy import logging as logg
from scipy.special import gammaln
from typing import Union, Dict, Optional

from dandelion.utilities._core import *
from dandelion.utilities._io import *
from dandelion.utilities._utilities import *
from dandelion.tools._chao1 import chao1
from dandelion.tools._gini import gini_index
from dandelion.tools._network import (
    clone_centrality,
    clone_degree,
    generate_network,
)
from dandelion.tools._shannon import shannon


def clone_rarefaction(
    self: Union[Dandelion, AnnData],
    groupby: str,
    clone_key: Optional[str] = None,
    diversity_key: Optional[str] = None,
    verbose: bool = False,
) -> Union[AnnData, Dict]:
    """
    Return rarefaction predictions for cell numbers vs clone size.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to split the calculation of clone numbers for a given number of cells for e.g. sample, patient etc.
    clone_key : str, Optional
        Column name specifying the clone_id column in metadata/obs.
    diversity_key : str, Optional
        key for 'diversity' results in AnnData's `.uns`.
    verbose: bool
        whether to print progess.

    Returns
    -------
    Dictionary containing rarefaction results or updated `.uns` slot if `AnnData` object is used.
    """
    start = logg.info("Constructing rarefaction curve")

    if isinstance(self, AnnData):
        metadata = self.obs.copy()
    elif isinstance(self, Dandelion):
        metadata = self.metadata.copy()

    if clone_key is None:
        clonekey = "clone_id"
    else:
        clonekey = clone_key

    groups = list(set(metadata[groupby]))
    metadata = metadata[metadata["contig_QC_pass"].isin([True, "True"])]
    if type(metadata[clonekey]) == "category":
        metadata[clonekey] = metadata[clonekey].cat.remove_unused_categories()
    res = {}
    for g in groups:
        _metadata = metadata[metadata[groupby] == g]
        res[g] = _metadata[clonekey].value_counts()
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
                np.array(
                    res_.iloc[
                        i,
                    ]
                ),
                z,
            )
            for z in n
        ]
    y = pd.DataFrame([rarecurve[c] for c in rarecurve]).T
    pred = pd.DataFrame(
        [np.append(np.arange(1, s, 10), s) for s in res_.sum(axis=1)],
        index=res_.index,
    ).T

    if diversity_key is None:
        diversitykey = "diversity"
    else:
        diversitykey = diversity_key

    if isinstance(self, AnnData):
        if diversitykey not in self.uns:
            self.uns[diversitykey] = {}
        self.uns[diversitykey] = {
            "rarefaction_cells_x": pred,
            "rarefaction_clones_y": y,
        }
    logg.info(
        " finished",
        time=start,
        deep=("updated `.uns` with rarefaction curves.\n"),
    )
    if isinstance(self, Dandelion):
        return {"rarefaction_cells_x": pred, "rarefaction_clones_y": y}


def clone_diversity(
    self: Union[Dandelion, AnnData],
    groupby: str,
    method: Literal["gini", "chao1", "shannon"] = "gini",
    metric: Literal["clone_network", "clone_degree", "clone_centrality"] = None,
    clone_key: Optional[str] = None,
    update_obs_meta: bool = True,
    diversity_key: Optional[str] = None,
    resample: bool = False,
    downsample: Optional[int] = None,
    n_resample: int = 50,
    normalize: bool = True,
    reconstruct_network: bool = True,
    expanded_only: bool = False,
    use_contracted: bool = False,
    key_added: Optional[str] = None,
    verbose: bool = False,
    **kwargs,
) -> Union[pd.DataFrame, Dandelion, AnnData]:
    """
    Compute B cell clones diversity : Gini indices, Chao1 estimates, or Shannon entropy.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the gini indices on, for e.g. sample, patient etc.
    method : str
        Method for diversity estimation. Either one of ['gini', 'chao1', 'shannon'].
    metric : str, Optional
        Metric to use for calculating Gini indices of clones.
        Accepts one of ['clone_network', 'clone_degree', 'clone_centrality'].
        `None` defaults to 'clone_network'.
    clone_key : str, Optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned.
        If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str, Optional
        key for 'diversity' results in `.uns`.
    downsample : int, Optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    resample : bool
        Whether or not to randomly sample cells without replacement to
        the minimum size of groups for the diversity calculation. Default is False.
    n_resample : int
        Number of times to perform resampling. Default is 50.
    normalize : bool
        Whether or not to return normalized Shannon Entropy according to https://math.stackexchange.com/a/945172.
        Default is True.
    reconstruct_network : bool
        Whether or not to reconstruct the network for Gini Index based measures.
        Default is True and will reconstruct for each group specified by groupby option.
    expanded_only : bool
        Whether or not to calculate gini indices using expanded clones only. Default is False i.e. use all cells/clones.
    use_contracted : bool
        Whether or not to perform the gini calculation after contraction of clone network.
        Only applies to calculation of clone size gini index. Default is False.
        This is to try and preserve the single-cell properties of the network.
    key_added : str, list, Optional
        column names for output.
    verbose : bool
        whether to print progress.
    **kwargs
        passed to dandelion.tl.generate_nework
    Returns
    -------
    `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    if downsample is not None:
        resample = True

    if method == "gini":
        if update_obs_meta:
            diversity_gini(
                self,
                groupby=groupby,
                metric=metric,
                clone_key=clone_key,
                update_obs_meta=update_obs_meta,
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
            return diversity_gini(
                self,
                groupby=groupby,
                metric=metric,
                clone_key=clone_key,
                update_obs_meta=update_obs_meta,
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
        if update_obs_meta:
            diversity_chao1(
                self,
                groupby=groupby,
                clone_key=clone_key,
                update_obs_meta=update_obs_meta,
                diversity_key=diversity_key,
                resample=resample,
                n_resample=n_resample,
                downsample=downsample,
                key_added=key_added,
                verbose=verbose,
            )
        else:
            return diversity_chao1(
                self,
                groupby=groupby,
                clone_key=clone_key,
                update_obs_meta=update_obs_meta,
                diversity_key=diversity_key,
                resample=resample,
                n_resample=n_resample,
                downsample=downsample,
                key_added=key_added,
                verbose=verbose,
            )
    if method == "shannon":
        if update_obs_meta:
            diversity_shannon(
                self,
                groupby=groupby,
                clone_key=clone_key,
                update_obs_meta=update_obs_meta,
                diversity_key=diversity_key,
                resample=resample,
                n_resample=n_resample,
                normalize=normalize,
                downsample=downsample,
                key_added=key_added,
                verbose=verbose,
            )
        else:
            return diversity_shannon(
                self,
                groupby=groupby,
                clone_key=clone_key,
                update_obs_meta=update_obs_meta,
                diversity_key=diversity_key,
                resample=resample,
                n_resample=n_resample,
                normalize=normalize,
                downsample=downsample,
                key_added=key_added,
                verbose=verbose,
            )


def clone_networkstats(
    self: Dandelion,
    expanded_only: bool = False,
    network_clustersize: bool = False,
    verbose: bool = False,
):
    """Retrieve network stats."""
    start = logg.info("Calculating vertex size of nodes after contraction")

    if isinstance(self, Dandelion):
        if self.graph is None:
            raise AttributeError(
                "Graph not found. Please run tl.generate_network."
            )
        else:
            if expanded_only:
                G = self.graph[1]
            else:
                G = self.graph[0]
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
    self: Union[Dandelion, AnnData],
    groupby: str,
    metric: Optional[str] = None,
    clone_key: Optional[str] = None,
    update_obs_meta: bool = False,
    diversity_key: Optional[str] = None,
    resample: bool = False,
    n_resample: int = 50,
    downsample: Optional[int] = None,
    reconstruct_network: bool = True,
    expanded_only: bool = False,
    use_contracted: bool = False,
    key_added: Optional[str] = None,
    verbose: bool = False,
    **kwargs,
) -> Union[pd.DataFrame, Dandelion]:
    """
    Compute B cell clones Gini indices.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Gini indices on, for e.g. sample, patient etc.
    metric : str, Optional
        Metric to use for calculating Gini indices of clones.
        Accepts one of ['clone_network', 'clone_degree', 'clone_centrality'].
        Defaults to 'clone_centrality'.
    clone_key : str, Optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned.
        If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str, Optional
        Key for 'diversity' results in `.uns`.
    resample : bool
        Whether or not to randomly sample cells without replacement to
        the minimum size of groups for the diversity calculation.
        Default is False. Resampling will automatically trigger reconstruction of network.
    n_resample : int
        Number of times to perform resampling. Default is 50.
    downsample : int, Optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    reconstruct_network : bool
        Whether or not to reconstruct the network for Gini Index based measures.
        Default is True and will reconstruct for each group specified by groupby option.
    expanded_only : bool
        Whether or not to calculate gini indices using expanded clones only. Default is False i.e. use all cells/clones.
    use_contracted : bool
        Whether or not to perform the gini calculation after contraction of clone network.
        Only applies to calculation of clone size gini index. Default is False.
        This is to try and preserve the single-cell properties of the network.
    key_added : str, list, Optional
        column names for output.
    verbose : bool
        whether to print progress.
    **kwargs
        passed to dandelion.tl.generate_nework
    Returns
    -------
    `pandas` dataframe or `Dandelion` object with updated `.metadata` slot.
    """
    start = logg.info("Calculating Gini indices")

    def gini_indices(
        self: Dandelion,
        groupby: str,
        metric: Optional[str] = None,
        clone_key: Optional[str] = None,
        resample: bool = False,
        n_resample: int = 50,
        downsample: Optional[int] = None,
        reconstruct_network: bool = True,
        expanded_only: bool = False,
        contracted: bool = False,
        key_added: Optional[str] = None,
        verbose: bool = False,
        **kwargs,
    ) -> pd.DataFrame:
        """Gini indices."""
        if isinstance(self, AnnData):
            raise TypeError("Only Dandelion class object accepted.")
        elif isinstance(self, Dandelion):
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = "clone_id"
        else:
            clonekey = clone_key

        if metric is None:
            met = "clone_network"
        else:
            met = metric

        # split up the table by groupby
        metadata[groupby] = metadata[groupby].astype("category")
        metadata[groupby].cat.remove_unused_categories(inplace=True)
        groups = list(set(metadata[groupby]))

        if downsample is None:
            minsize = metadata[groupby].value_counts().min()
        else:
            minsize = downsample
            if minsize > metadata[groupby].value_counts().min():
                logg.info(
                    "Downsampling size provided of {} was larger than the smallest group size. ".format(
                        downsample
                    )
                    + "Defaulting to the smallest group size for downsampling."
                )
                minsize = metadata[groupby].value_counts().min()

        if minsize < 100:
            logg.info(
                "The minimum cell numbers when grouped by {} is {}.".format(
                    groupby, minsize
                )
                + " Exercise caution when interpreting diversity measures."
            )

        res1 = {}
        if met == "clone_network":
            logg.info(
                "Computing Gini indices for cluster and vertex size using network."
            )
            if not reconstruct_network:
                n_n, v_s, c_s = clone_networkstats(
                    self,
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
                    g_c_v[vs] = gini_index(v_sizes, method="trapezoids")
                    if g_c_v[vs] < 0 or np.isnan(g_c_v[vs]):
                        g_c_v[vs] = 0
                    for cell in n_n:
                        g_c_v_res.update({cell: g_c_v[n_n[cell]]})
                c_sizes = np.array(
                    np.array(sorted(list(flatten(c_s.values())), reverse=True))
                )
                if len(c_sizes) > 1:
                    c_sizes = np.append(c_sizes, 0)
                g_c_c = gini_index(c_sizes, method="trapezoids")
                if g_c_c < 0 or np.isnan(g_c_c):
                    g_c_c = 0
                for cell in n_n:
                    g_c_c_res.update({cell: g_c_c})
                self.metadata["clone_network_vertex_size_gini"] = pd.Series(
                    g_c_v_res
                )
                self.metadata["clone_network_cluster_size_gini"] = pd.Series(
                    g_c_c_res
                )
        elif met == "clone_centrality":
            logg.info(
                "Computing gini indices for clone size using metadata and node closeness centrality using network."
            )
            clone_centrality(self, verbose=verbose)
        elif met == "clone_degree":
            logg.info(
                "Computing gini indices for clone size using metadata and node degree using network."
            )
            clone_degree(self, verbose=verbose)
        metadata = self.metadata.copy()
        data = self.data.copy()
        res2 = {}

        if resample:
            logg.info(
                "Downsampling each group specified in `{}` to {} cells for calculating gini indices.".format(
                    groupby, minsize
                )
            )

        for g in groups:
            # clone size distribution
            _dat = metadata[metadata[groupby] == g]
            _data = data[data["cell_id"].isin(list(_dat.index))]
            ddl_dat = Dandelion(_data, metadata=_dat)
            if resample:
                sizelist = []
                if isinstance(self, Dandelion):
                    graphlist = []
                for i in tqdm(
                    range(0, n_resample),
                    bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
                    disable=not verbose,
                ):
                    if isinstance(self, Dandelion):
                        resampled = generate_network(
                            ddl_dat,
                            clone_key=clonekey,
                            downsample=minsize,
                            verbose=verbose,
                            compute_layout=False,
                            **kwargs,
                        )
                        if met == "clone_network":
                            n_n, v_s, c_s = clone_networkstats(
                                resampled,
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
                                g_c_v[vs] = gini_index(
                                    v_sizes, method="trapezoids"
                                )
                                if g_c_v[vs] < 0 or np.isnan(g_c_v[vs]):
                                    g_c_v[vs] = 0
                                for cell in n_n:
                                    g_c_v_res.update({cell: g_c_v[n_n[cell]]})
                            c_sizes = np.array(
                                sorted(
                                    list(flatten(c_s.values())), reverse=True
                                )
                            )
                            if len(c_sizes) > 1:
                                c_sizes = np.append(c_sizes, 0)
                            g_c_c = gini_index(c_sizes, method="trapezoids")
                            if g_c_c < 0 or np.isnan(g_c_c):
                                g_c_c = 0
                            for cell in n_n:
                                g_c_c_res.update({cell: g_c_c})
                            resampled.metadata[
                                "clone_network_vertex_size_gini"
                            ] = pd.Series(g_c_v_res)
                            resampled.metadata[
                                "clone_network_cluster_size_gini"
                            ] = pd.Series(g_c_c_res)
                        elif met == "clone_centrality":
                            clone_centrality(resampled, verbose=verbose)
                        elif met == "clone_degree":
                            clone_degree(resampled, verbose=verbose)
                        else:
                            raise ValueError(
                                (
                                    "Unknown metric for calculating network stats. Please specify "
                                    + "one of `clone_network`, `clone_centrality` or `clone_degree`."
                                )
                            )
                        # clone size gini
                        _dat = resampled.metadata.copy()
                        _tab = _dat[clonekey].value_counts()
                        if "nan" in _tab.index or np.nan in _tab.index:
                            try:
                                _tab.drop("nan", inplace=True)
                            except:
                                _tab.drop(np.nan, inplace=True)
                        if met == "clone_network":
                            sizelist.append(
                                _dat[met + "_cluster_size_gini"].mean()
                            )
                        else:
                            clonesizecounts = np.array(_tab)
                            clonesizecounts = clonesizecounts[
                                clonesizecounts > 0
                            ]
                            if len(clonesizecounts) > 1:
                                # append a single zero for lorenz curve calculation
                                clonesizecounts = np.append(clonesizecounts, 0)
                            if len(clonesizecounts) > 0:
                                g_c = gini_index(
                                    clonesizecounts, method="trapezoids"
                                )
                                # probably not needed anymore but keep just in case
                                if g_c < 0 or np.isnan(g_c):
                                    g_c = 0
                            else:
                                g_c = 0
                            sizelist.append(g_c)

                        if met == "clone_network":
                            graphlist.append(
                                _dat[met + "_vertex_size_gini"].mean()
                            )
                        else:
                            # vertex closeness centrality or weighted degree distribution
                            # only calculate for expanded clones. If including non-expanded clones, the centrality is just zero which doesn't help.
                            connectednodes = resampled.metadata[met][
                                resampled.metadata[met] > 0
                            ]
                            graphcounts = np.array(
                                connectednodes.value_counts()
                            )
                            # graphcounts = np.append(graphcounts, 0) # if I add a  zero here, it will skew the results when the centrality measure is uniform.... so leave it out for now.
                            if len(graphcounts) > 0:
                                g_c = gini_index(
                                    graphcounts, method="trapezoids"
                                )
                                if g_c < 0 or np.isnan(g_c):
                                    g_c = 0
                            else:
                                g_c = 0
                            graphlist.append(g_c)
                try:
                    g_c = sum(sizelist) / len(sizelist)
                except:
                    g_c = 0

                res1.update({g: g_c})
                if "graphlist" in locals():
                    g_c = sum(graphlist) / len(graphlist)
                    try:
                        g_c = sum(graphlist) / len(graphlist)
                    except:
                        g_c = 0
                    res2.update({g: g_c})
            else:
                _tab = _dat[clonekey].value_counts()
                if "nan" in _tab.index or np.nan in _tab.index:
                    try:
                        _tab.drop("nan", inplace=True)
                    except:
                        _tab.drop(np.nan, inplace=True)
                if met != "clone_network":
                    clonesizecounts = np.array(_tab)
                    clonesizecounts = clonesizecounts[clonesizecounts > 0]
                    if len(clonesizecounts) > 1:
                        # append a single zero for lorenz curve calculation
                        clonesizecounts = np.append(clonesizecounts, 0)
                    if len(clonesizecounts) > 0:
                        g_c = gini_index(clonesizecounts, method="trapezoids")
                        # probably not needed anymore but keep just in case
                        if g_c < 0 or np.isnan(g_c):
                            g_c = 0
                    else:
                        g_c = 0
                    res1.update({g: g_c})
                if isinstance(self, Dandelion):
                    if met == "clone_network":
                        if reconstruct_network:
                            generate_network(
                                ddl_dat,
                                clone_key=clonekey,
                                verbose=verbose,
                                compute_layout=False,
                                **kwargs,
                            )
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
                                g_c_v[vs] = gini_index(
                                    v_sizes, method="trapezoids"
                                )
                                if g_c_v[vs] < 0 or np.isnan(g_c_v[vs]):
                                    g_c_v[vs] = 0
                                for cell in n_n:
                                    g_c_v_res.update({cell: g_c_v[n_n[cell]]})
                            c_sizes = np.array(
                                sorted(
                                    list(flatten(c_s.values())), reverse=True
                                )
                            )
                            if len(c_sizes) > 1:
                                c_sizes = np.append(c_sizes, 0)
                            g_c_c = gini_index(c_sizes, method="trapezoids")
                            if g_c_c < 0 or np.isnan(g_c_c):
                                g_c_c = 0
                            for cell in n_n:
                                g_c_c_res.update({cell: g_c_c})
                            # ddl_dat.metadata['clone_network_vertex_size_gini'] = pd.Series(g_c_v_res)
                            # ddl_dat.metadata['clone_network_cluster_size_gini'] = pd.Series(g_c_c_res)
                            res2.update({g: pd.Series(g_c_v_res).mean()})
                            res1.update({g: pd.Series(g_c_c_res).mean()})
                        else:
                            res2.update(
                                {g: _dat[met + "_vertex_size_gini"].mean()}
                            )
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
                        if len(graphcounts) > 0:
                            g_c = gini_index(graphcounts, method="trapezoids")
                            if g_c < 0 or np.isnan(g_c):
                                g_c = 0
                        else:
                            g_c = 0
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

    def transfer_gini_indices(
        self: Dandelion, gini_results: pd.DataFrame, groupby: str
    ) -> Dandelion:
        """Transfer gini indicies."""
        metadata = self.metadata.copy()

        groups = list(set(metadata[groupby]))
        for c in gini_results.columns:
            metadata[c] = np.nan
            for g in groups:
                for i in metadata.index:
                    if metadata.at[i, groupby] == g:
                        metadata.at[i, c] = gini_results[c][g]
        self.metadata = metadata.copy()

    res = gini_indices(
        self,
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

    if diversity_key is None:
        diversitykey = "diversity"
    else:
        diversitykey = diversity_key

    if update_obs_meta:
        res_ = res.copy()
        transfer_gini_indices(self, res_, groupby)
        if isinstance(self, Dandelion):
            logg.info(
                " finished",
                time=start,
                deep=("updated `.metadata` with Gini indices.\n"),
            )
    else:
        res_ = res.copy()
        logg.info(" finished", time=start)
        return res_


def diversity_chao1(
    self: Union[Dandelion, AnnData],
    groupby: str,
    clone_key: Optional[str] = None,
    update_obs_meta: bool = False,
    diversity_key: Optional[str] = None,
    resample: bool = False,
    n_resample: int = 50,
    downsample: Optional[int] = None,
    key_added: Optional[str] = None,
    verbose: bool = False,
) -> Union[pd.DataFrame, Dandelion, AnnData]:
    """
    Compute B cell clones Chao1 estimates.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Chao1 estimates on, for e.g. sample, patient etc.
    clone_key : str, Optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's
        metadata/obs slot.
    diversity_key : str, Optional
        key for 'diversity' results in `.uns`.
    resample : bool
        Whether or not to randomly sample cells without replacement to the minimum size of groups for the diversity
        calculation. Default is False.
    n_resample : int
        Number of times to perform resampling. Default is 50.
    downsample : int, Optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    key_added : str, list, Optional
        column names for output.
    verbose : bool
        whether to print progress.

    Returns
    -------
    `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    start = logg.info("Calculating Chao1 estimates")

    def chao1_estimates(
        self: Union[Dandelion, AnnData],
        groupby: str,
        clone_key: Optional[str] = None,
        resample: bool = False,
        n_resample: int = 50,
        downsample: Optional[int] = None,
        key_added: Optional[str] = None,
        verbose: bool = False,
    ) -> pd.DataFrame:
        """Chao1 estimates."""
        if isinstance(self, AnnData):
            metadata = self.obs.copy()
        elif isinstance(self, Dandelion):
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = "clone_id"
        else:
            clonekey = clone_key

        # split up the table by groupby
        metadata[groupby] = metadata[groupby].astype("category")
        metadata[groupby].cat.remove_unused_categories(inplace=True)
        groups = list(set(metadata[groupby]))

        if downsample is None:
            minsize = metadata[groupby].value_counts().min()
        else:
            minsize = downsample
            if minsize > metadata[groupby].value_counts().min():
                logg.info(
                    "Downsampling size provided of {} was larger than the smallest group size. ".format(
                        downsample
                    )
                    + "Defaulting to the smallest group size for downsampling."
                )
                minsize = metadata[groupby].value_counts().min()

        if minsize < 100:
            logg.info(
                "The minimum cell numbers when grouped by {} is {}.".format(
                    groupby, minsize
                )
                + " Exercise caution when interpreting diversity measures."
            )

        if resample:
            logg.info(
                "Downsampling each group specified in `{}` to {} cells for calculating Chao1 estimates.".format(
                    groupby, minsize
                )
            )
        res1 = {}
        for g in groups:
            # clone size distribution
            _dat = metadata[metadata[groupby] == g]
            if resample:
                sizelist = []
                graphlist = []
                for i in tqdm(
                    range(0, n_resample),
                    bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
                    disable=not verbose,
                ):
                    _dat = _dat.sample(minsize)
                    _tab = _dat[clonekey].value_counts()
                    if "nan" in _tab.index or np.nan in _tab.index:
                        try:
                            _tab.drop("nan", inplace=True)
                        except:
                            _tab.drop(np.nan, inplace=True)
                    clonesizecounts = np.array(_tab)
                    clonesizecounts = clonesizecounts[clonesizecounts > 0]
                    if len(clonesizecounts) > 0:
                        g_c = chao1(clonesizecounts)
                    else:
                        g_c = 0
                    sizelist.append(g_c)
                try:
                    g_c = sum(sizelist) / len(sizelist)
                except:
                    g_c = 0
                res1.update({g: g_c})
            else:
                _tab = _dat[clonekey].value_counts()
                if "nan" in _tab.index or np.nan in _tab.index:
                    try:
                        _tab.drop("nan", inplace=True)
                    except:
                        _tab.drop(np.nan, inplace=True)
                clonesizecounts = np.array(_tab)
                clonesizecounts = clonesizecounts[clonesizecounts > 0]
                if len(clonesizecounts) > 0:
                    g_c = chao1(clonesizecounts)
                else:
                    g_c = 0
                res1.update({g: g_c})

        res_df = pd.DataFrame.from_dict([res1]).T
        if key_added is None:
            res_df.columns = ["clone_size_chao1"]
        else:
            if type(key_added) is list:
                res_df.columns = key_added[0]
            else:
                res_df.columns = [key_added]

        return res_df

    def transfer_chao1_estimates(
        self: Union[Dandelion, AnnData],
        chao1_results: pd.DataFrame,
        groupby: str,
    ):
        """Transfer chao1 estimates."""
        if isinstance(self, AnnData):
            metadata = self.obs.copy()
        elif isinstance(self, Dandelion):
            metadata = self.metadata.copy()

        groups = list(set(metadata[groupby]))
        for c in chao1_results.columns:
            metadata[c] = np.nan
            for g in groups:
                for i in metadata.index:
                    if metadata.at[i, groupby] == g:
                        metadata.at[i, c] = chao1_results[c][g]
        if isinstance(self, AnnData):
            self.obs = metadata.copy()
        elif isinstance(self, Dandelion):
            self.metadata = metadata.copy()

    res = chao1_estimates(
        self,
        groupby=groupby,
        clone_key=clone_key,
        resample=resample,
        n_resample=n_resample,
        downsample=downsample,
        verbose=verbose,
    )

    if diversity_key is None:
        diversitykey = "diversity"
    else:
        diversitykey = diversity_key

    if isinstance(self, AnnData):
        if diversitykey not in self.uns:
            self.uns[diversitykey] = {}
        self.uns[diversitykey].update({"chao1": res})

    if update_obs_meta:
        res_ = res.copy()
        transfer_chao1_estimates(self, res_, groupby)
        if isinstance(self, Dandelion):
            logg.info(
                " finished",
                time=start,
                deep=("updated `.metadata` with Chao1 estimates.\n"),
            )
        elif isinstance(self, AnnData):
            logg.info(
                " finished",
                time=start,
                deep=("updated `.obs` and `.uns` with Chao1 estimates.\n"),
            )
    else:
        res_ = res.copy()
        if isinstance(self, AnnData):
            logg.info(
                " finished",
                time=start,
                deep=("updated `.uns` with Chao1 estimates.\n"),
            )
        else:
            logg.info(" finished", time=start)
        return res_


def diversity_shannon(
    self: Union[Dandelion, AnnData],
    groupby: str,
    clone_key: Optional[str] = None,
    update_obs_meta: bool = False,
    diversity_key: Optional[str] = None,
    resample: bool = False,
    n_resample: int = 50,
    normalize: bool = True,
    downsample: Optional[int] = None,
    key_added: Optional[str] = None,
    verbose: bool = False,
) -> Union[pd.DataFrame, Dandelion, AnnData]:
    """
    Compute B cell clones Shannon entropy.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Shannon entropy on, for e.g. sample, patient etc.
    clone_key : str, Optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str, Optional
        key for 'diversity' results in `.uns`.
    resample : bool
        Whether or not to randomly sample cells without replacement to the minimum size of groups for the diversity calculation. Default is False.
    n_resample : int
        Number of times to perform resampling. Default is 50.
    normalize : bool
        Whether or not to return normalized Shannon Entropy according to https://math.stackexchange.com/a/945172. Default is True.
    downsample : int, Optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    key_added : str, list, Optional
        column names for output.
    verbose : bool
        whether to print progress.
    Returns
    -------
    `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    start = logg.info("Calculating Shannon entropy")

    def shannon_entropy(
        self: Union[Dandelion, AnnData],
        groupby: str,
        clone_key: Optional[str] = None,
        resample: bool = False,
        n_resample: int = 50,
        normalize: bool = True,
        downsample: Optional[int] = None,
        key_added: Optional[str] = None,
        verbose: bool = False,
    ) -> pd.DataFrame:
        """Shannon entropy."""
        if isinstance(self, AnnData):
            metadata = self.obs.copy()
        elif isinstance(self, Dandelion):
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = "clone_id"
        else:
            clonekey = clone_key

        # split up the table by groupby
        metadata[groupby] = metadata[groupby].astype("category")
        metadata[groupby] = metadata[groupby].cat.remove_unused_categories()
        groups = list(set(metadata[groupby]))

        if downsample is None:
            minsize = metadata[groupby].value_counts().min()
        else:
            minsize = downsample
            if minsize > metadata[groupby].value_counts().min():
                logg.info(
                    "Downsampling size provided of {} was larger than the smallest group size. ".format(
                        downsample
                    )
                    + "Defaulting to the smallest group size for downsampling."
                )
                minsize = metadata[groupby].value_counts().min()

        if minsize < 100:
            logg.info(
                "The minimum cell numbers when grouped by {} is {}.".format(
                    groupby, minsize
                )
                + " Exercise caution when interpreting diversity measures."
            )

        if resample:
            logg.info(
                "Downsampling each group specified in `{}` to {} cells for calculating Shannon entropy.".format(
                    groupby, minsize
                )
            )

        res1 = {}
        sleep(0.5)
        for g in groups:
            # clone size distribution
            _dat = metadata[metadata[groupby] == g]
            if resample:
                sizelist = []
                graphlist = []

                for i in tqdm(
                    range(0, n_resample),
                    bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
                    disable=not verbose,
                ):
                    _dat = _dat.sample(minsize)
                    _tab = _dat[clonekey].value_counts()
                    if "nan" in _tab.index or np.nan in _tab.index:
                        try:
                            _tab.drop("nan", inplace=True)
                        except:
                            _tab.drop(np.nan, inplace=True)
                    clonesizecounts = np.array(_tab)
                    clonesizecounts = clonesizecounts[clonesizecounts > 0]
                    if len(clonesizecounts) > 0:
                        if normalize:
                            if len(clonesizecounts) == 1:
                                g_c = 0
                            else:
                                clonesizecounts_freqs = (
                                    clonesizecounts / np.sum(clonesizecounts)
                                )
                                g_c = -np.sum(
                                    (
                                        clonesizecounts_freqs
                                        * np.log(clonesizecounts_freqs)
                                    )
                                    / np.log(len(clonesizecounts_freqs))
                                )
                        else:
                            g_c = shannon(clonesizecounts)
                    else:
                        if normalize:
                            g_c = 1
                        else:
                            g_c = 0
                    sizelist.append(g_c)
                try:
                    g_c = sum(sizelist) / len(sizelist)
                except:
                    g_c = 0
                res1.update({g: g_c})
            else:
                _tab = _dat[clonekey].value_counts()
                if "nan" in _tab.index or np.nan in _tab.index:
                    try:
                        _tab.drop("nan", inplace=True)
                    except:
                        _tab.drop(np.nan, inplace=True)
                clonesizecounts = np.array(_tab)
                clonesizecounts = clonesizecounts[clonesizecounts > 0]
                if len(clonesizecounts) > 0:
                    if normalize:
                        if len(clonesizecounts) == 1:
                            g_c = 0
                        else:
                            clonesizecounts_freqs = clonesizecounts / np.sum(
                                clonesizecounts
                            )
                            g_c = -np.sum(
                                (
                                    clonesizecounts_freqs
                                    * np.log(clonesizecounts_freqs)
                                )
                                / np.log(len(clonesizecounts_freqs))
                            )
                    else:
                        g_c = shannon(clonesizecounts)
                else:
                    if normalize:
                        g_c = 1
                    else:
                        g_c = 0
                res1.update({g: g_c})

        res_df = pd.DataFrame.from_dict([res1]).T
        if key_added is None:
            if normalize:
                res_df.columns = ["clone_size_normalized_shannon"]
            else:
                res_df.columns = ["clone_size_shannon"]
        else:
            if type(key_added) is list:
                res_df.columns = key_added[0]
            else:
                res_df.columns = [key_added]

        return res_df

    def transfer_shannon_entropy(
        self: Union[Dandelion, AnnData],
        shannon_results: pd.DataFrame,
        groupby: str,
    ):
        """Transfer shannon entropy."""
        if isinstance(self, AnnData):
            metadata = self.obs.copy()
        elif isinstance(self, Dandelion):
            metadata = self.metadata.copy()

        groups = list(set(metadata[groupby]))
        for c in shannon_results.columns:
            metadata[c] = np.nan
            for g in groups:
                for i in metadata.index:
                    if metadata.at[i, groupby] == g:
                        metadata.at[i, c] = shannon_results[c][g]
        if isinstance(self, AnnData):
            self.obs = metadata.copy()
        elif isinstance(self, Dandelion):
            self.metadata = metadata.copy()

    res = shannon_entropy(
        self,
        groupby=groupby,
        clone_key=clone_key,
        resample=resample,
        n_resample=n_resample,
        normalize=normalize,
        downsample=downsample,
    )

    if diversity_key is None:
        diversitykey = "diversity"
    else:
        diversitykey = diversity_key

    if isinstance(self, AnnData):
        if diversitykey not in self.uns:
            self.uns[diversitykey] = {}
        self.uns[diversitykey].update({"shannon": res})

    if update_obs_meta:
        res_ = res.copy()
        transfer_shannon_entropy(self, res_, groupby)
        if isinstance(self, Dandelion):
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
        elif isinstance(self, AnnData):
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
    else:
        res_ = res.copy()
        if isinstance(self, AnnData):
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
        return res_


def chooseln(N, k):
    """
    R's lchoose in python
    from https://stackoverflow.com/questions/21767690/python-log-n-choose-k
    """
    return gammaln(N + 1) - gammaln(N - k + 1) - gammaln(k + 1)


def rarefun(y, sample):
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
