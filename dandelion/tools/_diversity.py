#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-08-13 21:08:53
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-08-14 16:26:37

import pandas as pd
import numpy as np
import networkx as nx
from ..utilities._utilities import *
from scipy.special import comb
from anndata import AnnData
from skbio.diversity.alpha import chao1, gini_index, shannon

def rarefaction(self, groupby, clone_key = None, diversity_key = None):
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
        key for 'diversity' results in `.uns`.
    Returns
    ----------
        `Dandelion` object with updated `.uns` slot or `AnnData` object with updated `.uns` slot with 'rarefaction` dictionary.
    """
    start = logg.info('Constructing rarefaction curves')

    if self.__class__ == AnnData:
        metadata = self.obs.copy()
    elif self.__class__ == Dandelion:
        metadata = self.metadata.copy()
    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    groups = list(set(metadata[groupby]))
    res = {}
    for g in groups:
        _metadata = metadata[metadata[groupby]==g]
        res[g] = _metadata[clonekey].value_counts()
    res_ = pd.DataFrame.from_dict(res, orient = 'index')
    res_.fillna(0, inplace = True)
    res_.reset_index(inplace = True)    
    
    pred, yhat, label = [], [], []
    _ = res_.apply(rarefactionCurve, args=(pred, yhat, label), axis=1)
    cells_pred = pd.DataFrame(pred, index = label)
    clones_yhat = pd.DataFrame(yhat, index = label)
    if diversity_key is None:
        diversitykey = 'diversity'
    else:
        diversitykey = diversity_key

    if diversitykey not in self.uns:
        self.uns[diversitykey] = {}
    self.uns[diversitykey] = {'rarefaction_cells_x':cells_pred, 'rarefaction_clones_y':clones_yhat}
    logg.info(' finished', time=start,
            deep=('updated `.uns` with rarefaction curves.\n'))

def rarefactionCurve(x, pred, yhat, label):
    # modified from ecopy: https://github.com/Auerilas/ecopy/blob/master/ecopy/diversity/rarefy.py
    ix = x['index']
    x_ = x.drop('index').astype('float')
    notabs = ~np.isnan(x_)
    y = x_[notabs]
    n = np.sum(y)
    Sn = len(x_)
    cell_Pred = np.linspace(0, n, 1000)
    clone_yHat = [rareCurve_Func(i, Sn, n, y) for i in cell_Pred]
    pred.append(cell_Pred)
    yhat.append(clone_yHat)
    label.append(str(ix))
    
def rareCurve_Func(i, Sn, n, x):
    # from ecopy: https://github.com/Auerilas/ecopy/blob/master/ecopy/diversity/rarefy.py
    sBar = Sn -  np.sum(comb(n-x, i))/comb(n, i) 
    return sBar


def clone_diversity(self, groupby, method = 'gini', splitby = None, clone_key = None, update_obs_meta = True, diversity_key = None):
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
    splitby : str, optional
        Column name to split by when calculating gini indices. None does not results in splitting.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str, optional
        key for 'diversity' results in `.uns`.
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    if method == 'gini':
        diversity_gini(self, groupby, splitby, clone_key, update_obs_meta, diversity_key)
    if method == 'chao1':
        diversity_chao1(self, groupby, splitby, clone_key, update_obs_meta, diversity_key)
    if method == 'shannon':
        diversity_shannon(self, groupby, splitby, clone_key, update_obs_meta, diversity_key)


def diversity_gini(self, groupby, splitby = None, clone_key = None, update_obs_meta = False, diversity_key = None):
    """
    Compute B cell clones Gini indices.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Gini indices on, for e.g. sample, patient etc.
    splitby : str, optional
        Column name to split by when calculating Gini indices. None does not results in splitting.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str, optional
        key for 'diversity' results in `.uns`.
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    start = logg.info('Calculating Gini indices')

    def gini_indices(self, groupby, splitby = None, clone_key = None):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = 'clone_id'
        else:
            clonekey = clone_key

        if 'clone_degree' not in metadata.columns:
            raise ValueError("`clone_degree` not found in provided object. Please run tl.clone_degree")
        # split up the table by groupby
        groups = list(set(metadata[groupby]))
        if splitby is None:
            res1 = {}
            res2 = {}
            for g in groups:
                # clone size distribution
                _dat = metadata[metadata[groupby] == g]
                clonesizecounts = np.array(_dat[clonekey].value_counts())
                if len(clonesizecounts) > 0:
                    g_c = gini_index(clonesizecounts)
                else:
                    g_c = np.nan
                res1.update({g:g_c})

                # vertex weighted degree distribution
                graphcounts = np.array(_dat['clone_degree'].value_counts())
                if len(graphcounts) > 0:
                    g_c = gini_index(graphcounts)
                else:
                    g_c = np.nan
                res2.update({g:g_c})
            res_df = pd.DataFrame.from_dict([res1,res2]).T
            res_df.columns = ['clone_size_gini', 'clone_degree_gini']
            return(res_df)
        else:
            splits = list(set(metadata[splitby]))
            res1 = Tree()
            res2 = Tree()
            for g in groups:
                for s in splits:
                    # clone size distribution
                    _dat = metadata[(metadata[groupby] == g) & (metadata[splitby] == s)]
                    clonesizecounts = np.array(_dat[clonekey].value_counts())
                    if len(clonesizecounts) > 0:
                        g_c = gini_index(clonesizecounts)
                    else:
                        g_c = np.nan
                    res1[g][s] = g_c

                    # vertex weighted degree distribution
                    graphcounts = np.array(_dat['clone_degree'].value_counts())
                    if len(graphcounts) > 0:
                        g_c = gini_index(graphcounts)
                    else:
                        g_c = np.nan
                    res2[g][s] = g_c
            res_df1 = pd.DataFrame([dict(res1[r]) for r in res1], index = res1.keys())
            res_df2 = pd.DataFrame([dict(res2[r]) for r in res2], index = res2.keys())
            res_df1.columns = ['clone_size_gini' + x for x in res_df1.columns]
            res_df2.columns = ['clone_degree_gini' + x for x in res_df2.columns]            
            return(res_df1,res_df2)

    def transfer_gini_indices(self, gini_results, groupby, splitby = None):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()

        groups = list(set(metadata[groupby]))
        if splitby is None:
            if not isinstance(gini_results, list):
                for c in gini_results.columns:
                    metadata[c] = np.nan
                    for g in groups:
                        for i in metadata.index:
                            if metadata.at[i, groupby] == g:
                                metadata.at[i, c] = gini_results[c][g]
        else:            
            splits = list(set(metadata[splitby]))
            for r in gini_results:
                for c in r.columns:
                    metadata[c] = np.nan
                    for g in groups:
                        for s in splits:
                            for i in metadata.index:
                                if metadata.at[i, groupby] == g:
                                    if metadata.at[i, splitby] == s:
                                        metadata.at[i, c] = r[c][g]
        if self.__class__ == AnnData:
            self.obs = metadata.copy()
        elif self.__class__ == Dandelion:
            self.metadata = metadata.copy()

    if diversity_key is None:
        diversitykey = 'diversity'
    else:
        diversitykey = diversity_key

    if diversitykey not in self.uns:
        self.uns[diversitykey] = {}
    if splitby is None:
        res  = gini_indices(self, groupby, splitby, clone_key)
    else:
        res1, res2  = gini_indices(self, groupby, splitby, clone_key)
        res = res1.join(res2)
    self.uns[diversitykey].update({'gini':res})

    if update_obs_meta:
        if splitby is None:
            res_ = res.copy()
        else:
            res_ = [res1, res2]
        transfer_gini_indices(self, res_, groupby, splitby)
        if self.__class__ == Dandelion:
            logg.info(' finished', time=start,
                deep=('updated `.metadata` and `.uns` with Gini indices.\n'))
        if self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.obs` and `.uns` with Gini indices.\n'))
    else:
        logg.info(' finished', time=start,
            deep=('updated `.uns` with Gini indices.\n'))

def diversity_chao1(self, groupby, splitby = None, clone_key = None, update_obs_meta = False, diversity_key = None):
    """
    Compute B cell clones Chao1 estimates.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Chao1 estimates on, for e.g. sample, patient etc.
    splitby : str, optional
        Column name to split by when calculating Chao1 estimates. None does not results in splitting.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str, optional
        key for 'diversity' results in `.uns`.
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    start = logg.info('Calculating Chao1 estimates')

    def chao1_estimates(self, groupby, splitby = None, clone_key = None):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = 'clone_id'
        else:
            clonekey = clone_key

        if 'clone_degree' not in metadata.columns:
            raise ValueError("`clone_degree` not found in provided object. Please run tl.clone_degree")
        # split up the table by groupby
        groups = list(set(metadata[groupby]))
        if splitby is None:
            res1 = {}
            res2 = {}
            for g in groups:
                # clone size distribution
                _dat = metadata[metadata[groupby] == g]
                clonesizecounts = np.array(_dat[clonekey].value_counts())
                if len(clonesizecounts) > 0:
                    g_c = chao1(clonesizecounts)
                else:
                    g_c = np.nan
                res1.update({g:g_c})

                # vertex weighted degree distribution
                graphcounts = np.array(_dat['clone_degree'].value_counts())
                if len(graphcounts) > 0:
                    g_c = chao1(graphcounts)
                else:
                    g_c = np.nan
                res2.update({g:g_c})
            res_df = pd.DataFrame.from_dict([res1,res2]).T
            res_df.columns = ['clone_size_chao1', 'clone_degree_chao1']
            return(res_df)
        else:
            splits = list(set(metadata[splitby]))
            res1 = Tree()
            res2 = Tree()
            for g in groups:
                for s in splits:
                    # clone size distribution
                    _dat = metadata[(metadata[groupby] == g) & (metadata[splitby] == s)]
                    clonesizecounts = np.array(_dat[clonekey].value_counts())
                    if len(clonesizecounts) > 0:
                        g_c = chao1(clonesizecounts)
                    else:
                        g_c = np.nan
                    res1[g][s] = g_c

                    # vertex weighted degree distribution
                    graphcounts = np.array(_dat['clone_degree'].value_counts())
                    if len(graphcounts) > 0:
                        g_c = chao1(graphcounts)
                    else:
                        g_c = np.nan
                    res2[g][s] = g_c
            res_df1 = pd.DataFrame([dict(res1[r]) for r in res1], index = res1.keys())
            res_df2 = pd.DataFrame([dict(res2[r]) for r in res2], index = res2.keys())
            res_df1.columns = ['clone_size_chao1' + x for x in res_df1.columns]
            res_df2.columns = ['clone_degree_chao1' + x for x in res_df2.columns]            
            return(res_df1,res_df2)

    def transfer_chao1_estimates(self, chao1_results, groupby, splitby = None):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()

        groups = list(set(metadata[groupby]))
        if splitby is None:
            if not isinstance(chao1_results, list):
                for c in chao1_results.columns:
                    metadata[c] = np.nan
                    for g in groups:
                        for i in metadata.index:
                            if metadata.at[i, groupby] == g:
                                metadata.at[i, c] = chao1_results[c][g]
        else:            
            splits = list(set(metadata[splitby]))
            for r in chao1_results:
                for c in r.columns:
                    metadata[c] = np.nan
                    for g in groups:
                        for s in splits:
                            for i in metadata.index:
                                if metadata.at[i, groupby] == g:
                                    if metadata.at[i, splitby] == s:
                                        metadata.at[i, c] = r[c][g]
        if self.__class__ == AnnData:
            self.obs = metadata.copy()
        elif self.__class__ == Dandelion:
            self.metadata = metadata.copy()

    if diversity_key is None:
        diversitykey = 'diversity'
    else:
        diversitykey = diversity_key

    if diversitykey not in self.uns:
        self.uns[diversitykey] = {}
    if splitby is None:
        res  = chao1_estimates(self, groupby, splitby, clone_key)
    else:
        res1, res2  = chao1_estimates(self, groupby, splitby, clone_key)
        res = res1.join(res2)
    self.uns[diversitykey].update({'chao1':res})

    if update_obs_meta:
        if splitby is None:
            res_ = res.copy()
        else:
            res_ = [res1, res2]
        transfer_chao1_estimates(self, res_, groupby, splitby)
        if self.__class__ == Dandelion:
            logg.info(' finished', time=start,
                deep=('updated `.metadata` and `.uns` with Chao1 estimates.\n'))
        if self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.obs` and `.uns` with Chao1 estimates.\n'))
    else:
        logg.info(' finished', time=start,
            deep=('updated `.uns` with Chao1 estimates.\n'))

def diversity_shannon(self, groupby, splitby = None, clone_key = None, update_obs_meta = False, diversity_key = None):
    """
    Compute B cell clones Shannon entropy.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        Column name to calculate the Shannon entropy on, for e.g. sample, patient etc.
    splitby : str, optional
        Column name to split by when calculating Shannon entropy. None does not results in splitting.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata.
    update_obs_meta : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
    diversity_key : str, optional
        key for 'diversity' results in `.uns`.
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    start = logg.info('Calculating Shannon entropy')

    def shannon_entropy(self, groupby, splitby = None, clone_key = None):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()
        if clone_key is None:
            clonekey = 'clone_id'
        else:
            clonekey = clone_key

        if 'clone_degree' not in metadata.columns:
            raise ValueError("`clone_degree` not found in provided object. Please run tl.clone_degree")
        # split up the table by groupby
        groups = list(set(metadata[groupby]))
        if splitby is None:
            res1 = {}
            res2 = {}
            for g in groups:
                # clone size distribution
                _dat = metadata[metadata[groupby] == g]
                clonesizecounts = np.array(_dat[clonekey].value_counts())
                if len(clonesizecounts) > 0:
                    g_c = shannon(clonesizecounts)
                else:
                    g_c = np.nan
                res1.update({g:g_c})

                # vertex weighted degree distribution
                graphcounts = np.array(_dat['clone_degree'].value_counts())
                if len(graphcounts) > 0:
                    g_c = shannon(graphcounts)
                else:
                    g_c = np.nan
                res2.update({g:g_c})
            res_df = pd.DataFrame.from_dict([res1,res2]).T
            res_df.columns = ['clone_size_shannon', 'clone_degree_shannon']
            return(res_df)
        else:
            splits = list(set(metadata[splitby]))
            res1 = Tree()
            res2 = Tree()
            for g in groups:
                for s in splits:
                    # clone size distribution
                    _dat = metadata[(metadata[groupby] == g) & (metadata[splitby] == s)]
                    clonesizecounts = np.array(_dat[clonekey].value_counts())
                    if len(clonesizecounts) > 0:
                        g_c = shannon(clonesizecounts)
                    else:
                        g_c = np.nan
                    res1[g][s] = g_c

                    # vertex weighted degree distribution
                    graphcounts = np.array(_dat['clone_degree'].value_counts())
                    if len(graphcounts) > 0:
                        g_c = shannon(graphcounts)
                    else:
                        g_c = np.nan
                    res2[g][s] = g_c
            res_df1 = pd.DataFrame([dict(res1[r]) for r in res1], index = res1.keys())
            res_df2 = pd.DataFrame([dict(res2[r]) for r in res2], index = res2.keys())
            res_df1.columns = ['clone_size_shannon' + x for x in res_df1.columns]
            res_df2.columns = ['clone_degree_shannon' + x for x in res_df2.columns]            
            return(res_df1,res_df2)

    def transfer_shannon_entropy(self, shannon_results, groupby, splitby = None):
        if self.__class__ == AnnData:
            metadata = self.obs.copy()
        elif self.__class__ == Dandelion:
            metadata = self.metadata.copy()

        groups = list(set(metadata[groupby]))
        if splitby is None:
            if not isinstance(shannon_results, list):
                for c in shannon_results.columns:
                    metadata[c] = np.nan
                    for g in groups:
                        for i in metadata.index:
                            if metadata.at[i, groupby] == g:
                                metadata.at[i, c] = shannon_results[c][g]
        else:            
            splits = list(set(metadata[splitby]))
            for r in shannon_results:
                for c in r.columns:
                    metadata[c] = np.nan
                    for g in groups:
                        for s in splits:
                            for i in metadata.index:
                                if metadata.at[i, groupby] == g:
                                    if metadata.at[i, splitby] == s:
                                        metadata.at[i, c] = r[c][g]
        if self.__class__ == AnnData:
            self.obs = metadata.copy()
        elif self.__class__ == Dandelion:
            self.metadata = metadata.copy()

    if diversity_key is None:
        diversitykey = 'diversity'
    else:
        diversitykey = diversity_key

    if diversitykey not in self.uns:
        self.uns[diversitykey] = {}
    if splitby is None:
        res  = shannon_entropy(self, groupby, splitby, clone_key)
    else:
        res1, res2  = shannon_entropy(self, groupby, splitby, clone_key)
        res = res1.join(res2)
    self.uns[diversitykey].update({'shannon':res})

    if update_obs_meta:
        if splitby is None:
            res_ = res.copy()
        else:
            res_ = [res1, res2]
        transfer_shannon_entropy(self, res_, groupby, splitby)
        if self.__class__ == Dandelion:
            logg.info(' finished', time=start,
                deep=('updated `.metadata` and `.uns` with Shannon entropy.\n'))
        if self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.obs` and `.uns` with Shannon entropy.\n'))
    else:
        logg.info(' finished', time=start,
            deep=('updated `.uns` with Shannon entropy.\n'))
