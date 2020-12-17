#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-08-13 21:08:53
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-12-17 14:02:13

import pandas as pd
import numpy as np
import networkx as nx
from ..utilities._utilities import *
from ..tools._network import clone_centrality, clone_degree, generate_network
from scipy.special import gammaln
from anndata import AnnData
from skbio.diversity.alpha import chao1, gini_index, shannon
from tqdm import tqdm
from time import sleep
import warnings

def clone_rarefaction(self, groupby, clone_key=None, diversity_key = None):
    """
    Returns rarefaction predictions for cell numbers vs clone size.
    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to split the calculation of clone numbers for a given number of cells for e.g. sample, patient etc.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata/obs.
    diversity_key : str, optional
        key for 'diversity' results in AnnData's `.uns`.
    Returns
    ----------
        Dictionary containing rarefaction results or updated `.uns` slot if `AnnData` object is used.
    """
    start = logg.info('Constructing rarefaction curve')

    if self.__class__ == AnnData:
        metadata = self.obs.copy()
    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    groups = list(set(metadata[groupby]))
    metadata = metadata[metadata['bcr_QC_pass'].isin([True, 'True'])]
    metadata[clonekey] = metadata[clonekey].cat.remove_unused_categories()
    res = {}
    for g in groups:
        _metadata = metadata[metadata[groupby]==g]
        res[g] = _metadata[clonekey].value_counts()
    res_ = pd.DataFrame.from_dict(res, orient = 'index')

    # remove those with no counts
    rowsum = res_.sum(axis = 1)
    print('removing due to zero counts:', ', '.join([res_.index[i] for i, x in enumerate(res_.sum(axis = 1) == 0) if x]))
    res_ = res_[~(res_.sum(axis = 1) == 0)]

    # set up for calculating rarefaction
    tot = res_.apply(sum, axis = 1)
    S = res_.apply(lambda x: x[x > 0].shape[0], axis = 1)
    nr = res_.shape[0]

    # append the results to a dictionary
    rarecurve = {}
    sleep(0.5)
    for i in tqdm(range(0, nr), desc = 'Calculating rarefaction curve '):
        n = np.arange(1, tot[i], step = 10)
        if (n[-1:] != tot[i]):
            n = np.append(n, tot[i])
        rarecurve[res_.index[i]] = [rarefun(np.array(res_.iloc[i,]), z) for z in n]
    y = pd.DataFrame([rarecurve[c] for c in rarecurve]).T
    pred = pd.DataFrame([np.append(np.arange(1, s, 10),s) for s in res_.sum(axis = 1)], index = res_.index).T

    if diversity_key is None:
        diversitykey = 'diversity'
    else:
        diversitykey = diversity_key

    if self.__class__ == AnnData:
        if diversitykey not in self.uns:
            self.uns[diversitykey] = {}
        self.uns[diversitykey] = {'rarefaction_cells_x':pred, 'rarefaction_clones_y':y}
    logg.info(' finished', time=start,
            deep=('updated `.uns` with rarefaction curves.\n'))
    if self.__class__ == Dandelion:
        return({'rarefaction_cells_x':pred, 'rarefaction_clones_y':y})

def clone_diversity(self, groupby, method = 'gini', metric = None, clone_key = None, update_obs_meta = True, diversity_key = None, resample = False, downsample = None, n_resample = 50, normalize = True):
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
    metric : str, optional
        Metric to use for calculating Gini indices of clones. Accepts one of ['clone_vertexsize', clone_degree', 'clone_centrality']. Defaults to 'clone_vertexsize'.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str, optional
        key for 'diversity' results in `.uns`.
    downsample : int, optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    resample : bool
        Whether or not to randomly sample cells without replacement to the minimum size of groups for the diversity calculation. Default is False.
    n_resample : int
        Number of times to perform resampling. Default is 50.
    normalize : bool
        Whether or not to return normalized Shannon Entropy according to https://math.stackexchange.com/a/945172. Default is True.
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    if downsample is not None:
        resample = True
    if method == 'gini':
        if update_obs_meta:
            diversity_gini(self, groupby=groupby, metric=metric, clone_key=clone_key, update_obs_meta=update_obs_meta, diversity_key=diversity_key, resample=resample, n_resample=n_resample, downsample=downsample)
        else:
            return(diversity_gini(self, groupby=groupby, metric=metric, clone_key=clone_key, update_obs_meta=update_obs_meta, diversity_key=diversity_key, resample=resample, n_resample=n_resample, downsample=downsample))
    if method == 'chao1':
        if update_obs_meta:
            diversity_chao1(self, groupby=groupby, clone_key=clone_key, update_obs_meta=update_obs_meta, diversity_key=diversity_key, resample=resample, n_resample=n_resample, downsample=downsample)
        else:
            return(diversity_chao1(self, groupby=groupby, clone_key=clone_key, update_obs_meta=update_obs_meta, diversity_key=diversity_key, resample=resample, n_resample=n_resample, downsample=downsample))
    if method == 'shannon':
        if update_obs_meta:
            diversity_shannon(self, groupby=groupby, clone_key=clone_key, update_obs_meta=update_obs_meta, diversity_key=diversity_key, resample=resample, n_resample=n_resample, normalize = normalize, downsample=downsample)
        else:
            return(diversity_shannon(self, groupby=groupby, clone_key=clone_key, update_obs_meta=update_obs_meta, diversity_key=diversity_key, resample=resample, n_resample=n_resample, normalize = normalize, downsample=downsample))

def clone_vertexsize(self, verbose = True):
    if verbose:
        start = logg.info('Calculating vertex size of nodes after contraction')

    if self.__class__ == Dandelion:
        try:
            G = self.graph[0]
        except:
            dist = np.sum([self.distance[x].toarray() for x in self.distance if type(self.distance[x]) is csr_matrix], axis = 0)
            A = csr_matrix(dist)
            G = nx.Graph()
            G.add_weighted_edges_from(zip(list(self.metadata.index), list(self.metadata.index), A.data))

        if len(G) is 0:
            raise AttributeError('Graph not found. Plase run tl.generate_network.')
        else:
            remove_edges = defaultdict(list)
            vertexsizes = defaultdict(list)
            nodes_names = defaultdict(list)
            vertex_counts = defaultdict(dict)
            if verbose:
                for subg in tqdm(nx.connected_components(G), desc = 'Reducing graph '):
                    nodes = sorted(list(subg))
                    tmp = nodes[0] #  just assign the value in a single cell, because this will be representative of the clone
                    for n in nodes:
                        nodes_names[n] = tmp # keep so i can reference later
                    if len(nodes) > 1:
                        G_ = G.subgraph(nodes).copy()
                        remove_edges[tmp] = [(e[0],e[1]) for e in G_.edges(data = True) if e[2]['weight'] > 0]
                        if len(remove_edges[tmp]) > 0:
                            G_.remove_edges_from(remove_edges[tmp])
                            for connected in nx.connected_components(G_):
                                vertexsizes[tmp].append(len(connected))
                            vertexsizes[tmp] = sorted(vertexsizes[tmp], reverse = True)
                        else:
                            vertexsizes[tmp] = [1 for i in range(len(G_.edges(data = True)))]
                    else:
                        vertexsizes[tmp] = [1]
            else:
                for subg in nx.connected_components(G):
                    nodes = sorted(list(subg))
                    tmp = nodes[0] # just assign the value in a single cell, because this will be representative of the clone
                    for n in nodes:
                        nodes_names[n] = tmp # keep so i can reference later
                    if len(nodes) > 1:
                        G_ = G.subgraph(nodes).copy()
                        remove_edges[tmp] = [(e[0],e[1]) for e in G_.edges(data = True) if e[2]['weight'] > 0]
                        if len(remove_edges[tmp]) > 0:
                            G_.remove_edges_from(remove_edges[tmp])
                            for connected in nx.connected_components(G_):
                                vertexsizes[tmp].append(len(connected))
                            vertexsizes[tmp] = sorted(vertexsizes[tmp], reverse = True)
                        else:
                            vertexsizes[tmp] = [1 for i in range(len(G_.edges(data = True)))]
                    else:
                        vertexsizes[tmp] = [1]
            # vertexsizes = pd.DataFrame(pd.Series(vertex_counts).reset_index().set_axis(['cell_id','counts'],1,inplace=False))
            return(nodes_names, vertexsizes)
    else:
        raise TypeError('Input object must be of {}'.format(Dandelion))


def diversity_gini(self, groupby, metric = None, clone_key = None, update_obs_meta = False, diversity_key = None, resample = False, n_resample = 50, downsample = None):
    """
    Compute B cell clones Gini indices.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Gini indices on, for e.g. sample, patient etc.
    metric : str, optional
        Metric to use for calculating Gini indices of clones. Accepts 'clone_degree' and 'clone_centrality'. Defaults to 'clone_centrality'.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str, optional
        Key for 'diversity' results in `.uns`.
    resample : bool
        Whether or not to randomly sample cells without replacement to the minimum size of groups for the diversity calculation. Default is False.
    n_resample : int
        Number of times to perform resampling. Default is 50.
    downsample : int, optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    start = logg.info('Calculating Gini indices')

    def gini_indices(self, groupby, metric = None, clone_key = None, resample = False, n_resample = 50, downsample = None):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = 'clone_id'
        else:
            clonekey = clone_key

        if metric is None:
            met = 'clone_vertexsize'
        else:
            met = metric

        # split up the table by groupby
        metadata[groupby] = metadata[groupby].astype('category')
        metadata[groupby].cat.remove_unused_categories(inplace = True)
        groups = list(set(metadata[groupby]))

        if downsample is None:
            minsize = metadata[groupby].value_counts().min()
        else:
            minsize = downsample
            if minsize > metadata[groupby].value_counts().min():
                print('Downsampling size provided of {} was larger than the smallest group size. Defaulting to the smallest group size for downsampling.'.format(downsample))
                minsize = metadata[groupby].value_counts().min()

        if minsize < 100:
            warnings.warn("The minimum cell numbers when grouped by {} is {} (group {}). Exercise caution when interpreting diversity measures.".format(groupby, minsize, metadata[groupby].value_counts().idxmin()))

        res1 = {}
        if self.__class__ == Dandelion:
            print("{} provided. Computing Gini indices for clone size and clone network.".format(self.__class__.__name__))
            sleep(0.5)
            if met == 'clone_vertexsize':
                n_n, v_s = clone_vertexsize(self, verbose = True)
                g_c = defaultdict(dict)
                g_c_res = {}
                for vs in v_s:
                    v_sizes = np.array(v_s[vs])
                    if len(v_sizes) > 1:
                        v_sizes = np.append(v_sizes, 0)
                    g_c[vs] = gini_index(v_sizes, method = 'trapezoids')
                    if g_c[vs] < 0 or np.isnan(g_c[vs]):
                        g_c[vs] = 0
                    for cell in n_n:
                        g_c_res.update({cell:g_c[n_n[cell]]})
                self.metadata['clone_vertexsize_gini'] = pd.Series(g_c_res)
            elif met == 'clone_centrality':
                clone_centrality(self, verbose = True)
            elif met == 'clone_degree':
                clone_degree(self, verbose = True)
            metadata = self.metadata.copy()
            data = self.data.copy()
            res2 = {}
        else:
            print("{} provided. Only computing Gini indices for clone size.".format(self.__class__.__name__))
        if resample:
            print("Downsampling each group specified in `{}` to {} cells for calculating gini indices.".format(groupby, minsize))
        sleep(0.5)
        for g in groups:
            # clone size distribution
            _dat = metadata[metadata[groupby] == g]
            if self.__class__ == Dandelion:
                _data = data[data['cell_id'].isin(list(_dat.index))]
                ddl_dat = Dandelion(_data, metadata = _dat)
            if resample:
                sizelist = []
                if self.__class__ == Dandelion:
                    graphlist = []
                for i in tqdm(range(0, n_resample)):
                    if self.__class__ == Dandelion:
                        resampled = generate_network(ddl_dat, clone_key = clonekey, downsample = minsize, verbose = False)
                        if met == 'clone_vertexsize':
                            n_n, v_s = clone_vertexsize(resampled, verbose = False)
                            g_c = defaultdict(dict)
                            g_c_res = {}
                            for vs in v_s:
                                v_sizes = np.array(v_s[vs])
                                if len(v_sizes) > 1:
                                    v_sizes = np.append(v_sizes, 0)
                                g_c[vs] = gini_index(v_sizes, method = 'trapezoids')
                                if g_c[vs] < 0 or np.isnan(g_c[vs]):
                                    g_c[vs] = 0
                                for cell in n_n:
                                    g_c_res.update({cell:g_c[n_n[cell]]})
                            resampled.metadata['clone_vertexsize_gini'] = pd.Series(g_c_res)
                        elif met == 'clone_centrality':
                            clone_centrality(resampled, verbose = False)
                        elif met == 'clone_degree':
                            clone_degree(resampled, verbose = False)
                        else:
                            raise ValueError('Unknown metric for calculating network stats. Please specify one of `clone_centrality` or `clone_degree`.')
                        # clone size gini
                        _dat = resampled.metadata.copy()
                        _tab = _dat[clonekey].value_counts()
                        if 'nan' in _tab.index or np.nan in _tab.index:
                            try:
                                _tab.drop('nan', inplace = True)
                            except:
                                _tab.drop(np.nan, inplace = True)
                        clonesizecounts = np.array(_tab)
                        clonesizecounts = clonesizecounts[clonesizecounts > 0]
                        if len(clonesizecounts) > 1:
                            # append a single zero for lorenz curve calculation
                            clonesizecounts = np.append(clonesizecounts, 0)
                        if len(clonesizecounts) > 0:
                            g_c = gini_index(clonesizecounts, method = 'trapezoids')
                            if g_c < 0 or np.isnan(g_c): # probably not needed anymore but keep just in case
                                g_c = 0
                        else:
                            g_c = 0
                        sizelist.append(g_c)

                        if met == 'clone_vertexsize':
                            graphlist.append(_dat[met+'_gini'].mean())
                        else:
                            # vertex closeness centrality or weighted degree distribution
                            connectednodes = resampled.metadata[met][resampled.metadata[met] > 0] # only calculate for expanded clones. If including non-expanded clones, the centrality is just zero which doesn't help.
                            graphcounts = np.array(connectednodes.value_counts())
                            # graphcounts = np.append(graphcounts, 0) # if I add a  zero here, it will skew the results when the centrality measure is uniform.... so leave it out for now.
                            if len(graphcounts) > 0:
                                g_c = gini_index(graphcounts, method = 'trapezoids')
                                if g_c < 0 or np.isnan(g_c):
                                    g_c = 0
                            else:
                                g_c = 0
                            graphlist.append(g_c)
                    elif self.__class__ == AnnData:
                        _dat = _dat.sample(minsize)
                        _tab = _dat[clonekey].value_counts()
                        if 'nan' in _tab.index or np.nan in _tab.index:
                            try:
                                _tab.drop('nan', inplace = True)
                            except:
                                _tab.drop(np.nan, inplace = True)
                        clonesizecounts = np.array(_tab)
                        clonesizecounts = clonesizecounts[clonesizecounts > 0]
                        # append a single zero for lorenz curve calculation
                        if len(clonesizecounts) > 1:
                            # append a single zero for lorenz curve calculation
                            clonesizecounts = np.append(clonesizecounts, 0)
                        if len(clonesizecounts) > 0:
                            g_c = gini_index(clonesizecounts, method = 'trapezoids')
                            if g_c < 0 or np.isnan(g_c): # probably not needed anymore but keep just in case
                                g_c = 0
                        else:
                            g_c = 0
                        sizelist.append(g_c)
                try:
                    g_c = sum(sizelist)/len(sizelist)
                except:
                    g_c = 0

                res1.update({g:g_c})
                if 'graphlist' in locals():
                    g_c = sum(graphlist)/len(graphlist)
                    try:
                        g_c = sum(graphlist)/len(graphlist)
                    except:
                        g_c = 0
                    res2.update({g:g_c})
            else:
                _tab = _dat[clonekey].value_counts()
                if 'nan' in _tab.index or np.nan in _tab.index:
                    try:
                        _tab.drop('nan', inplace = True)
                    except:
                        _tab.drop(np.nan, inplace = True)
                clonesizecounts = np.array(_tab)
                clonesizecounts = clonesizecounts[clonesizecounts > 0]
                if len(clonesizecounts) > 1:
                    # append a single zero for lorenz curve calculation
                    clonesizecounts = np.append(clonesizecounts, 0)
                if len(clonesizecounts) > 0:
                    g_c = gini_index(clonesizecounts, method = 'trapezoids')
                    if g_c < 0 or np.isnan(g_c): # probably not needed anymore but keep just in case
                        g_c = 0
                else:
                    g_c = 0
                res1.update({g:g_c})

                if self.__class__ == Dandelion:
                    if met == 'clone_vertexsize':
                        res2.update({g:_dat[met+'_gini'].mean()})
                    else:
                        # vertex closeness centrality or weighted degree distribution
                        connectednodes = _dat[met][_dat[met] > 0] # only calculate for expanded clones. If including non-expanded clones, the centrality is just zero which doesn't help.
                        graphcounts = np.array(connectednodes.value_counts())
                        # graphcounts = np.append(graphcounts, 0) # if I add a  zero here, it will skew the results when the centrality measure is uniform.... so leave it out for now.
                        if len(graphcounts) > 0:
                            g_c = gini_index(graphcounts, method = 'trapezoids')
                            if g_c < 0 or np.isnan(g_c):
                                g_c = 0
                        else:
                            g_c = 0
                        res2.update({g:g_c})

        if 'res2' in locals():
            res_df = pd.DataFrame.from_dict([res1,res2]).T
            res_df.columns = ['clone_size_gini', met + '_gini']
        else:
            res_df = pd.DataFrame.from_dict([res1]).T
            res_df.columns = ['clone_size_gini']
        return(res_df)

    def transfer_gini_indices(self, gini_results, groupby):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()

        groups = list(set(metadata[groupby]))
        for c in gini_results.columns:
            metadata[c] = np.nan
            for g in groups:
                for i in metadata.index:
                    if metadata.at[i, groupby] == g:
                        metadata.at[i, c] = gini_results[c][g]
        if self.__class__ == AnnData:
            self.obs = metadata.copy()
        elif self.__class__ == Dandelion:
            self.metadata = metadata.copy()

    res  = gini_indices(self, groupby, clone_key, resample = resample, n_resample = n_resample, downsample = downsample)

    if diversity_key is None:
        diversitykey = 'diversity'
    else:
        diversitykey = diversity_key

    if self.__class__ == AnnData:
        if diversitykey not in self.uns:
            self.uns[diversitykey] = {}
        self.uns[diversitykey].update({'gini':res})

    if update_obs_meta:
        res_ = res.copy()
        transfer_gini_indices(self, res_, groupby)
        sleep(0.5)
        if self.__class__ == Dandelion:
            logg.info(' finished', time=start,
                deep=('updated `.metadata` with Gini indices.\n'))
        elif self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.obs` and `.uns` with Gini indices.\n'))
    else:
        res_ = res.copy()
        sleep(0.5)
        if self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.uns` with Gini indices.\n'))
        else:
            logg.info(' finished', time=start)
        return(res_)

def diversity_chao1(self, groupby, clone_key = None, update_obs_meta = False, diversity_key = None, resample = False, n_resample = 50, downsample = None):
    """
    Compute B cell clones Chao1 estimates.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Chao1 estimates on, for e.g. sample, patient etc.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str, optional
        key for 'diversity' results in `.uns`.
    resample : bool
        Whether or not to randomly sample cells without replacement to the minimum size of groups for the diversity calculation. Default is False.
    n_resample : int
        Number of times to perform resampling. Default is 50.
    downsample : int, optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    start = logg.info('Calculating Chao1 estimates')

    def chao1_estimates(self, groupby, clone_key = None, resample = False, n_resample = 50, downsample = None):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = 'clone_id'
        else:
            clonekey = clone_key

        # split up the table by groupby
        metadata[groupby] = metadata[groupby].astype('category')
        metadata[groupby].cat.remove_unused_categories(inplace = True)
        groups = list(set(metadata[groupby]))

        if downsample is None:
            minsize = metadata[groupby].value_counts().min()
        else:
            minsize = downsample
            if minsize > metadata[groupby].value_counts().min():
                print('Downsampling size provided of {} was larger than the smallest group size. Defaulting to the smallest group size for downsampling.'.format(downsample))
                minsize = metadata[groupby].value_counts().min()

        if minsize < 100:
            warnings.warn('The minimum cell numbers when grouped by {} is {}. Exercise caution when interpreting diversity measures.'.format(groupby, minsize))

        if resample:
            print("Downsampling each group specified in `{}` to {} cells for calculating Chao1 estimates.".format(groupby, minsize))
        res1 = {}
        sleep(0.5)
        for g in groups:
            # clone size distribution
            _dat = metadata[metadata[groupby] == g]
            if resample:
                sizelist = []
                graphlist = []
                for i in tqdm(range(0, n_resample)):
                    _dat = _dat.sample(minsize)
                    _tab = _dat[clonekey].value_counts()
                    if 'nan' in _tab.index or np.nan in _tab.index:
                        try:
                            _tab.drop('nan', inplace = True)
                        except:
                            _tab.drop(np.nan, inplace = True)
                    clonesizecounts = np.array(_tab)
                    clonesizecounts = clonesizecounts[clonesizecounts > 0]
                    if len(clonesizecounts) > 0:
                        g_c = chao1(clonesizecounts)
                    else:
                        g_c = 0
                    sizelist.append(g_c)
                try:
                    g_c = sum(sizelist)/len(sizelist)
                except:
                    g_c = 0
                res1.update({g:g_c})
            else:
                _tab = _dat[clonekey].value_counts()
                if 'nan' in _tab.index or np.nan in _tab.index:
                    try:
                        _tab.drop('nan', inplace = True)
                    except:
                        _tab.drop(np.nan, inplace = True)
                clonesizecounts = np.array(_tab)
                clonesizecounts = clonesizecounts[clonesizecounts > 0]
                if len(clonesizecounts) > 0:
                    g_c = chao1(clonesizecounts)
                else:
                    g_c = 0
                res1.update({g:g_c})

        res_df = pd.DataFrame.from_dict([res1]).T
        res_df.columns = ['clone_size_chao1']

        return(res_df)

    def transfer_chao1_estimates(self, chao1_results, groupby):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()

        groups = list(set(metadata[groupby]))
        for c in chao1_results.columns:
            metadata[c] = np.nan
            for g in groups:
                for i in metadata.index:
                    if metadata.at[i, groupby] == g:
                        metadata.at[i, c] = chao1_results[c][g]
        if self.__class__ == AnnData:
            self.obs = metadata.copy()
        elif self.__class__ == Dandelion:
            self.metadata = metadata.copy()

    res  = chao1_estimates(self, groupby, clone_key, resample = resample, n_resample = n_resample, downsample = downsample)

    if diversity_key is None:
        diversitykey = 'diversity'
    else:
        diversitykey = diversity_key

    if self.__class__ == AnnData:
        if diversitykey not in self.uns:
            self.uns[diversitykey] = {}
        self.uns[diversitykey].update({'chao1':res})

    if update_obs_meta:
        res_ = res.copy()
        transfer_chao1_estimates(self, res_, groupby)
        sleep(0.5)
        if self.__class__ == Dandelion:
            logg.info(' finished', time=start,
                deep=('updated `.metadata` with Chao1 estimates.\n'))
        elif self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.obs` and `.uns` with Chao1 estimates.\n'))
    else:
        res_ = res.copy()
        sleep(0.5)
        if self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.uns` with Chao1 estimates.\n'))
        else:
            logg.info(' finished', time=start)
        return(res_)

def diversity_shannon(self, groupby, clone_key = None, update_obs_meta = False, diversity_key = None, resample = False, n_resample = 50, normalize = True, downsample = None):
    """
    Compute B cell clones Shannon entropy.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Shannon entropy on, for e.g. sample, patient etc.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str, optional
        key for 'diversity' results in `.uns`.
    resample : bool
        Whether or not to randomly sample cells without replacement to the minimum size of groups for the diversity calculation. Default is False.
    n_resample : int
        Number of times to perform resampling. Default is 50.
    normalize : bool
        Whether or not to return normalized Shannon Entropy according to https://math.stackexchange.com/a/945172. Default is True.
    downsample : int, optional
        number of cells to downsample to. If None, defaults to size of smallest group.
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    start = logg.info('Calculating Shannon entropy')

    def shannon_entropy(self, groupby, clone_key = None, resample = False, n_resample = 50, normalize = True, downsample = None):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = 'clone_id'
        else:
            clonekey = clone_key

        # split up the table by groupby
        metadata[groupby] = metadata[groupby].astype('category')
        metadata[groupby].cat.remove_unused_categories(inplace = True)
        groups = list(set(metadata[groupby]))

        if downsample is None:
            minsize = metadata[groupby].value_counts().min()
        else:
            minsize = downsample
            if minsize > metadata[groupby].value_counts().min():
                print('Downsampling size provided of {} was larger than the smallest group size. Defaulting to the smallest group size for downsampling.'.format(downsample))
                minsize = metadata[groupby].value_counts().min()

        if minsize < 100:
            warnings.warn('The minimum cell numbers when grouped by {} is {}. Exercise caution when interpreting diversity measures.'.format(groupby, minsize))

        if resample:
            print("Downsampling each group specified in `{}` to {} cells for calculating Shannon entropy.".format(groupby, minsize))

        res1 = {}
        sleep(0.5)
        for g in groups:
            # clone size distribution
            _dat = metadata[metadata[groupby] == g]
            if resample:
                sizelist = []
                graphlist = []

                for i in tqdm(range(0, n_resample)):
                    _dat = _dat.sample(minsize)
                    _tab = _dat[clonekey].value_counts()
                    if 'nan' in _tab.index or np.nan in _tab.index:
                        try:
                            _tab.drop('nan', inplace = True)
                        except:
                            _tab.drop(np.nan, inplace = True)
                    clonesizecounts = np.array(_tab)
                    clonesizecounts = clonesizecounts[clonesizecounts > 0]
                    if len(clonesizecounts) > 0:
                        if normalize:
                            if len(clonesizecounts) == 1:
                                g_c = 0
                            else:
                                clonesizecounts_freqs = clonesizecounts / np.sum(clonesizecounts)
                                g_c = -np.sum((clonesizecounts_freqs * np.log(clonesizecounts_freqs)) / np.log(len(clonesizecounts_freqs)))
                        else:
                            g_c = shannon(clonesizecounts)
                    else:
                        if normalize:
                            g_c = 1
                        else:
                            g_c = 0
                    sizelist.append(g_c)
                try:
                    g_c = sum(sizelist)/len(sizelist)
                except:
                    g_c = 0
                res1.update({g:g_c})
            else:
                _tab = _dat[clonekey].value_counts()
                if 'nan' in _tab.index or np.nan in _tab.index:
                    try:
                        _tab.drop('nan', inplace = True)
                    except:
                        _tab.drop(np.nan, inplace = True)
                clonesizecounts = np.array(_tab)
                clonesizecounts = clonesizecounts[clonesizecounts > 0]
                if len(clonesizecounts) > 0:
                    if normalize:
                        if len(clonesizecounts) == 1:
                            g_c = 0
                        else:
                            clonesizecounts_freqs = clonesizecounts / np.sum(clonesizecounts)
                            g_c = -np.sum((clonesizecounts_freqs * np.log(clonesizecounts_freqs)) / np.log(len(clonesizecounts_freqs)))
                    else:
                        g_c = shannon(clonesizecounts)
                else:
                    if normalize:
                        g_c = 1
                    else:
                        g_c = 0
                res1.update({g:g_c})

        res_df = pd.DataFrame.from_dict([res1]).T
        if normalize:
            res_df.columns = ['clone_size_normalized_shannon']
        else:
            res_df.columns = ['clone_size_shannon']
        return(res_df)

    def transfer_shannon_entropy(self, shannon_results, groupby):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()

        groups = list(set(metadata[groupby]))
        for c in shannon_results.columns:
            metadata[c] = np.nan
            for g in groups:
                for i in metadata.index:
                    if metadata.at[i, groupby] == g:
                        metadata.at[i, c] = shannon_results[c][g]
        if self.__class__ == AnnData:
            self.obs = metadata.copy()
        elif self.__class__ == Dandelion:
            self.metadata = metadata.copy()

    res  = shannon_entropy(self, groupby, clone_key, resample = resample, n_resample = n_resample, normalize = normalize, downsample = downsample)

    if diversity_key is None:
        diversitykey = 'diversity'
    else:
        diversitykey = diversity_key

    if self.__class__ == AnnData:
        if diversitykey not in self.uns:
            self.uns[diversitykey] = {}
        self.uns[diversitykey].update({'shannon':res})

    if update_obs_meta:
        res_ = res.copy()
        transfer_shannon_entropy(self, res_, groupby)
        sleep(0.5)
        if self.__class__ == Dandelion:
            if normalize:
                logg.info(' finished', time=start,
                    deep=('updated `.metadata` with normalized Shannon entropy.\n'))
            else:
                logg.info(' finished', time=start,
                    deep=('updated `.metadata` with Shannon entropy.\n'))
        elif self.__class__ == AnnData:
            if normalize:
                logg.info(' finished', time=start,
                    deep=('updated `.obs` and `.uns` with normalized Shannon entropy.\n'))
            else:
                logg.info(' finished', time=start,
                    deep=('updated `.obs` and `.uns` with Shannon entropy.\n'))
    else:
        res_ = res.copy()
        sleep(0.5)
        if self.__class__ == AnnData:
            if normalize:
                logg.info(' finished', time=start,
                    deep=('updated `.uns` with normalized Shannon entropy.\n'))
            else:
                logg.info(' finished', time=start,
                    deep=('updated `.uns` with Shannon entropy.\n'))
        else:
            logg.info(' finished', time=start)
        return(res_)

def chooseln(N, k):
    '''
    R's lchoose in python
    from https://stackoverflow.com/questions/21767690/python-log-n-choose-k
    '''
    return gammaln(N+1) - gammaln(N-k+1) - gammaln(k+1)

def rarefun(y, sample):
    '''
    Adapted from rarefun from vegan:
    https://github.com/vegandevs/vegan/blob/master/R/rarefy.R
    '''
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
    return(out)