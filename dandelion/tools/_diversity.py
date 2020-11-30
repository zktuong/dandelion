#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-08-13 21:08:53
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-11-30 12:11:25

import pandas as pd
import numpy as np
import networkx as nx
from ..utilities._utilities import *
from scipy.special import gammaln
from anndata import AnnData
from skbio.diversity.alpha import chao1, gini_index, shannon
from tqdm import tqdm
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
    metadata = metadata[metadata['has_bcr'].isin([True, 'True'])]
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
    for i in range(0, nr):
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

def clone_diversity(self, groupby, method = 'gini', metric = None, clone_key = None, update_obs_meta = True, diversity_key = None, resample = False, n_resample = 50, normalize = True):
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
        Metric to use for calculating Gini indices of clones. Accepts 'clone_degree' and 'clone_centrality'. Defaults to 'clone_centrality'.    
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
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    if method == 'gini':
        if update_obs_meta:
            diversity_gini(self, groupby, metric, clone_key, update_obs_meta, diversity_key, resample, n_resample)
        else:
            return(diversity_gini(self, groupby, metric, clone_key, update_obs_meta, diversity_key, resample, n_resample))
    if method == 'chao1':
        if update_obs_meta:
            diversity_chao1(self, groupby, clone_key, update_obs_meta, diversity_key, resample, n_resample)
        else:
            return(diversity_chao1(self, groupby, clone_key, update_obs_meta, diversity_key, resample, n_resample))
    if method == 'shannon':
        if update_obs_meta:
            diversity_shannon(self, groupby, clone_key, update_obs_meta, diversity_key, resample, n_resample, normalize)
        else:
            return(diversity_shannon(self, groupby, clone_key, update_obs_meta, diversity_key, resample, n_resample, normalize))

def diversity_gini(self, groupby, metric = None, clone_key = None, update_obs_meta = False, diversity_key = None, resample = False, n_resample = 50):
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
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    start = logg.info('Calculating Gini indices')

    def gini_indices(self, groupby, metric = None, clone_key = None, resample = False, n_resample = 50):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = 'clone_id'
        else:
            clonekey = clone_key

        if metric is None:
            met = 'clone_centrality'
        else:
            met = metric
        
        if met not in metadata.columns:
            from ..tools._network import clone_centrality, clone_degree
            if met == 'clone_centrality':                
                raise ValueError("`clone_centrality` not found in metadata. Please run tl.clone_centrality")
            elif met == 'clone_degree':
                raise ValueError("`clone_degree` not found in metadata. Please run tl.clone_degree")
            else:
                raise ValueError("`metric` not recognised. Please specify either `clone_centrality` or `clone_degree`.")
        # split up the table by groupby
        groups = list(set(metadata[groupby]))
            
        minsize = metadata[groupby].value_counts().min()
        if minsize < 100:
            warnings.warn("The minimum cell numbers when grouped by {} is {}. Practise caution when interpreting diversity measures.".format(groupby, minsize))

        res1 = {}
        res2 = {}
        for g in groups:
            # clone size distribution
            _dat = metadata[metadata[groupby] == g]
            if resample:                
                sizelist = []
                graphlist = []
                for i in range(0, n_resample):
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
                        g_c = gini_index(clonesizecounts, method = 'trapezoids')
                    else:
                        g_c = np.nan
                    sizelist.append(g_c)
                
                    # vertex closeness centrality or weighted degree distribution
                    graphcounts = np.array(_dat[met].value_counts())
                    if len(graphcounts) > 0:
                        g_c = gini_index(graphcounts)
                    else:
                        g_c = np.nan
                    graphlist.append(g_c)                
                try:
                    g_c = sum(sizelist)/len(sizelist)
                except:
                    g_c = np.nan
                res1.update({g:g_c})
                g_c = sum(graphlist)/len(graphlist)
                try:
                    g_c = sum(graphlist)/len(graphlist)
                except:
                    g_c = np.nan
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
                if len(clonesizecounts) > 0:
                    g_c = gini_index(clonesizecounts)
                else:
                    g_c = np.nan
                res1.update({g:g_c})
                # vertex closeness centrality or weighted degree distribution
                graphcounts = np.array(_dat[met].value_counts())
                if len(graphcounts) > 0:
                    g_c = gini_index(graphcounts, method = 'trapezoids')
                else:
                    g_c = np.nan
                res2.update({g:g_c})

        res_df = pd.DataFrame.from_dict([res1,res2]).T
        res_df.columns = ['clone_size_gini', met + '_gini']
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

    res  = gini_indices(self, groupby, clone_key, resample = resample, n_resample = n_resample)

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
        if self.__class__ == Dandelion:
            logg.info(' finished', time=start,
                deep=('updated `.metadata` with Gini indices.\n'))
        elif self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.obs` and `.uns` with Gini indices.\n'))
    else:
        res_ = res.copy()                
        if self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.uns` with Gini indices.\n'))
        else:
            logg.info(' finished', time=start)
        return(res_)

def diversity_chao1(self, groupby, clone_key = None, update_obs_meta = False, diversity_key = None, resample = False, n_resample = 50):
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
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    start = logg.info('Calculating Chao1 estimates')

    def chao1_estimates(self, groupby, clone_key = None, resample = False, n_resample = 50):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = 'clone_id'
        else:
            clonekey = clone_key

        # split up the table by groupby
        groups = list(set(metadata[groupby]))
            
        minsize = metadata[groupby].value_counts().min()
        if minsize < 100:
            warnings.warn('The minimum cell numbers when grouped by {} is {}. Practise caution when interpreting diversity measures.'.format(groupby, minsize))

        res1 = {}
        res2 = {}
        for g in groups:
            # clone size distribution
            _dat = metadata[metadata[groupby] == g]
            if resample:                
                sizelist = []
                graphlist = []
                for i in range(0, n_resample):
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
                        g_c = np.nan
                    sizelist.append(g_c)
                try:
                    g_c = sum(sizelist)/len(sizelist)
                except:
                    g_c = np.nan
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
                    g_c = np.nan
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
    
    res  = chao1_estimates(self, groupby, clone_key, resample = resample, n_resample = n_resample)
    
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
        if self.__class__ == Dandelion:
            logg.info(' finished', time=start,
                deep=('updated `.metadata` with Chao1 estimates.\n'))
        elif self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.obs` and `.uns` with Chao1 estimates.\n'))
    else:
        res_ = res.copy()                
        if self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.uns` with Chao1 estimates.\n'))
        else:
            logg.info(' finished', time=start)
        return(res_)

def diversity_shannon(self, groupby, clone_key = None, update_obs_meta = False, diversity_key = None, resample = False, n_resample = 50, normalize = True):
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
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    start = logg.info('Calculating Shannon entropy')

    def shannon_entropy(self, groupby, clone_key = None, resample = False, n_resample = 50, normalize = True):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = 'clone_id'
        else:
            clonekey = clone_key

        # split up the table by groupby
        groups = list(set(metadata[groupby]))
            
        minsize = metadata[groupby].value_counts().min()
        if minsize < 100:
            warnings.warn('The minimum cell numbers when grouped by {} is {}. Practise caution when interpreting diversity measures.'.format(groupby, minsize))

        res1 = {}
        res2 = {}
        for g in groups:
            # clone size distribution
            _dat = metadata[metadata[groupby] == g]
            if resample:                
                sizelist = []
                graphlist = []
                for i in range(0, n_resample):
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
                        g_c = np.nan
                    sizelist.append(g_c)              
                try:
                    g_c = sum(sizelist)/len(sizelist)
                except:
                    g_c = np.nan
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
                    g_c = np.nan
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

    res  = shannon_entropy(self, groupby, clone_key, resample = resample, n_resample = n_resample, normalize = normalize)

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