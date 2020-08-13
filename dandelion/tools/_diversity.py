#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-08-13 21:08:53
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-08-13 21:38:27

import pandas as pd
import numpy as np
import networkx as nx
from ..utilities._utilities import *
from anndata import AnnData
from skbio.diversity.alpha import chao1, gini_index, shannon

def diversity(self, groupby, method = 'gini', splitby = None, clone_key = None, return_dataframe = True):
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
    return_dataframe : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
    Returns
    ----------
        `pandas` dataframe, `Dandelion` object with updated `.metadata` slot or `AnnData` object with updated `.obs` slot.
    """
    if method == 'gini':
        diversity_gini(self, groupby, splitby, clone_key, return_dataframe)
    if method == 'chao1':
        diversity_chao1(self, groupby, splitby, clone_key, return_dataframe)
    if method == 'shannon':
        diversity_shannon(self, groupby, splitby, clone_key, return_dataframe)


def diversity_gini(self, groupby, splitby = None, clone_key = None, return_dataframe = True):
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
    return_dataframe : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
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

        if 'degree' not in metadata:
            raise ValueError("`degree` not found in provided object. Please run tl.clone_degree")
        # split up the table by groupby
        groups = list(set(metadata[groupby]))
        if splitby is None:
            res1 = {}
            res2 = {}
            for g in groups:
                # clone size distribution
                _dat = metadata[metadata[groupby] == g]
                clonesizecounts = np.array(_dat[clonekey].value_counts())
                g_c = gini_index(clonesizecounts)
                res1.update({g:g_c})

                # vertex weighted degree distribution
                graphcounts = np.array(_dat['clone_degree'].value_counts())
                g_c = gini_index(graphcounts)
                res2.update({g:g_c})
            res_df = pd.DataFrame.from_dict([res1,res2]).T
            res_df.columns = ['clone_size_gini', 'clone_degree_gini']
        else:
            splits = list(set(metadata[splitby]))
            res1 = Tree()
            res2 = Tree()
            for g in groups:
                for s in splits:
                    # clone size distribution
                    _dat = metadata[(metadata[groupby] == g) & (metadata[splitby] == s)]
                    clonesizecounts = np.array(_dat[clonekey].value_counts())
                    g_c = gini_index(clonesizecounts)
                    res1[g][s] = g_c

                    # vertex weighted degree distribution
                    graphcounts = np.array(_dat['clone_degree'].value_counts())
                    g_c = gini_index(graphcounts)
                    res2[g][s] = g_c
            res_df1 = pd.DataFrame([dict(res1[r]) for r in res1], index = res1.keys())
            res_df2 = pd.DataFrame([dict(res2[r]) for r in res2], index = res2.keys())
            res_df1.columns = ['clone_size_gini' + x for x in res_df1.columns]
            res_df2.columns = ['clone_degree_gini' + x for x in res_df2.columns]
            res_df = res_df1.join(res_df2)
        return(res_df)

    def _gini_indices(self, groupby, splitby = None, clone_key = None):
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
        if splitby is None:
            res1 = {}
            res2 = {}
            for g in groups:
                # clone size distribution
                _dat = metadata[metadata[groupby] == g]
                clonesizecounts = np.array(_dat[clonekey].value_counts())
                g_c = gini_index(clonesizecounts)
                res1.update({g:g_c})

                # vertex weighted degree distribution
                graphcounts = np.array(_dat['clone_degree'].value_counts())
                g_c = gini_index(graphcounts)
                res2.update({g:g_c})
            res_df = pd.DataFrame.from_dict([res1,res2]).T
            res_df.columns = ['clone_size_gini', 'clone_degree_gini']
        else:
            splits = list(set(metadata[splitby]))
            res1 = Tree()
            res2 = Tree()
            for g in groups:
                for s in splits:
                    # clone size distribution
                    _dat = metadata[(metadata[groupby] == g) & (metadata[splitby] == s)]
                    clonesizecounts = np.array(_dat[clonekey].value_counts())
                    g_c = gini_index(clonesizecounts)
                    res1[g][s] = g_c

                    # vertex weighted degree distribution
                    graphcounts = np.array(_dat['clone_degree'].value_counts())
                    g_c = gini_index(graphcounts)
                    res2[g][s] = g_c
            res_df1 = pd.DataFrame([dict(res1[r]) for r in res1], index = res1.keys())
            res_df2 = pd.DataFrame([dict(res2[r]) for r in res2], index = res2.keys())
            res_df1.columns = ['clone_size_gini_' + x for x in res_df1.columns]
            res_df2.columns = ['clone_degree_gini_' + x for x in res_df2.columns]
            res_df = [res_df1, res_df2]
        return(res_df)

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

    if return_dataframe:
        return(gini_indices(self, groupby, splitby, clone_key))
        logg.info(' finished', time=start)
    else:
        results = _gini_indices(self, groupby, splitby, clone_key)
        transfer_gini_indices(self, results, groupby, splitby)
        if self.__class__ == Dandelion:
            logg.info(' finished', time=start,
                deep=('updated `.metadata` with Gini indices\n'))
        if self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.obs` with Gini indices\n'))

def diversity_chao1(self, groupby, splitby = None, clone_key = None, return_dataframe = True):
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
    return_dataframe : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
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

        if 'degree' not in metadata:
            raise ValueError("`degree` not found in provided object. Please run tl.clone_degree")
        # split up the table by groupby
        groups = list(set(metadata[groupby]))
        if splitby is None:
            res1 = {}
            res2 = {}
            for g in groups:
                # clone size distribution
                _dat = metadata[metadata[groupby] == g]
                clonesizecounts = np.array(_dat[clonekey].value_counts())
                g_c = chao1(clonesizecounts)
                res1.update({g:g_c})

                # vertex weighted degree distribution
                graphcounts = np.array(_dat['clone_degree'].value_counts())
                g_c = chao1(graphcounts)
                res2.update({g:g_c})
            res_df = pd.DataFrame.from_dict([res1,res2]).T
            res_df.columns = ['clone_size_chao1', 'clone_degree_chao1']
        else:
            splits = list(set(metadata[splitby]))
            res1 = Tree()
            res2 = Tree()
            for g in groups:
                for s in splits:
                    # clone size distribution
                    _dat = metadata[(metadata[groupby] == g) & (metadata[splitby] == s)]
                    clonesizecounts = np.array(_dat[clonekey].value_counts())
                    g_c = chao1(clonesizecounts)
                    res1[g][s] = g_c

                    # vertex weighted degree distribution
                    graphcounts = np.array(_dat['clone_degree'].value_counts())
                    g_c = chao1(graphcounts)
                    res2[g][s] = g_c
            res_df1 = pd.DataFrame([dict(res1[r]) for r in res1], index = res1.keys())
            res_df2 = pd.DataFrame([dict(res2[r]) for r in res2], index = res2.keys())
            res_df1.columns = ['clone_size_chao1' + x for x in res_df1.columns]
            res_df2.columns = ['clone_degree_chao1' + x for x in res_df2.columns]
            res_df = res_df1.join(res_df2)
        return(res_df)

    def _chao1(self, groupby, splitby = None, clone_key = None):
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
        if splitby is None:
            res1 = {}
            res2 = {}
            for g in groups:
                # clone size distribution
                _dat = metadata[metadata[groupby] == g]
                clonesizecounts = np.array(_dat[clonekey].value_counts())
                g_c = chao1(clonesizecounts)
                res1.update({g:g_c})

                # vertex weighted degree distribution
                graphcounts = np.array(_dat['clone_degree'].value_counts())
                g_c = chao1(graphcounts)
                res2.update({g:g_c})
            res_df = pd.DataFrame.from_dict([res1,res2]).T
            res_df.columns = ['clone_size_chao1', 'clone_degree_chao1']
        else:
            splits = list(set(metadata[splitby]))
            res1 = Tree()
            res2 = Tree()
            for g in groups:
                for s in splits:
                    # clone size distribution
                    _dat = metadata[(metadata[groupby] == g) & (metadata[splitby] == s)]
                    clonesizecounts = np.array(_dat[clonekey].value_counts())
                    g_c = chao1(clonesizecounts)
                    res1[g][s] = g_c

                    # vertex weighted degree distribution
                    graphcounts = np.array(_dat['clone_degree'].value_counts())
                    g_c = chao1(graphcounts)
                    res2[g][s] = g_c
            res_df1 = pd.DataFrame([dict(res1[r]) for r in res1], index = res1.keys())
            res_df2 = pd.DataFrame([dict(res2[r]) for r in res2], index = res2.keys())
            res_df1.columns = ['clone_size_chao1_' + x for x in res_df1.columns]
            res_df2.columns = ['clone_degree_chao1_' + x for x in res_df2.columns]
            res_df = [res_df1, res_df2]
        return(res_df)

    def transfer_chao1(self, chao1_results, groupby, splitby = None):
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

    if return_dataframe:
        return(chao1_estimates(self, groupby, splitby, clone_key))
        logg.info(' finished', time=start)
    else:
        results = _chao1_estimates(self, groupby, splitby, clone_key)
        transfer_chao1(self, results, groupby, splitby)
        if self.__class__ == Dandelion:
            logg.info(' finished', time=start,
                deep=('updated `.metadata` with Chao1 estimates\n'))
        if self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.obs` with Chao1 estimates\n'))

def diversity_shannon(self, groupby, splitby = None, clone_key = None, return_dataframe = True):
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
    return_dataframe : bool
        If True, a `pandas` dataframe is returned. If False, function will try to populate the input object's metadata/obs slot.
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

        if 'degree' not in metadata:
            raise ValueError("`degree` not found in provided object. Please run tl.clone_degree")
        # split up the table by groupby
        groups = list(set(metadata[groupby]))
        if splitby is None:
            res1 = {}
            res2 = {}
            for g in groups:
                # clone size distribution
                _dat = metadata[metadata[groupby] == g]
                clonesizecounts = np.array(_dat[clonekey].value_counts())
                g_c = shannon(clonesizecounts)
                res1.update({g:g_c})

                # vertex weighted degree distribution
                graphcounts = np.array(_dat['clone_degree'].value_counts())
                g_c = shannon(graphcounts)
                res2.update({g:g_c})
            res_df = pd.DataFrame.from_dict([res1,res2]).T
            res_df.columns = ['clone_size_shannon', 'clone_degree_shannon']
        else:
            splits = list(set(metadata[splitby]))
            res1 = Tree()
            res2 = Tree()
            for g in groups:
                for s in splits:
                    # clone size distribution
                    _dat = metadata[(metadata[groupby] == g) & (metadata[splitby] == s)]
                    clonesizecounts = np.array(_dat[clonekey].value_counts())
                    g_c = shannon(clonesizecounts)
                    res1[g][s] = g_c

                    # vertex weighted degree distribution
                    graphcounts = np.array(_dat['clone_degree'].value_counts())
                    g_c = shannon(graphcounts)
                    res2[g][s] = g_c
            res_df1 = pd.DataFrame([dict(res1[r]) for r in res1], index = res1.keys())
            res_df2 = pd.DataFrame([dict(res2[r]) for r in res2], index = res2.keys())
            res_df1.columns = ['clone_size_shannon' + x for x in res_df1.columns]
            res_df2.columns = ['clone_degree_shannon' + x for x in res_df2.columns]
            res_df = res_df1.join(res_df2)
        return(res_df)

    def _shannon_entropy(self, groupby, splitby = None, clone_key = None):
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
        if splitby is None:
            res1 = {}
            res2 = {}
            for g in groups:
                # clone size distribution
                _dat = metadata[metadata[groupby] == g]
                clonesizecounts = np.array(_dat[clonekey].value_counts())
                g_c = shannon(clonesizecounts)
                res1.update({g:g_c})

                # vertex weighted degree distribution
                graphcounts = np.array(_dat['clone_degree'].value_counts())
                g_c = shannon(graphcounts)
                res2.update({g:g_c})
            res_df = pd.DataFrame.from_dict([res1,res2]).T
            res_df.columns = ['clone_size_shannon', 'clone_degree_shannon']
        else:
            splits = list(set(metadata[splitby]))
            res1 = Tree()
            res2 = Tree()
            for g in groups:
                for s in splits:
                    # clone size distribution
                    _dat = metadata[(metadata[groupby] == g) & (metadata[splitby] == s)]
                    clonesizecounts = np.array(_dat[clonekey].value_counts())
                    g_c = shannon(clonesizecounts)
                    res1[g][s] = g_c

                    # vertex weighted degree distribution
                    graphcounts = np.array(_dat['clone_degree'].value_counts())
                    g_c = shannon(graphcounts)
                    res2[g][s] = g_c
            res_df1 = pd.DataFrame([dict(res1[r]) for r in res1], index = res1.keys())
            res_df2 = pd.DataFrame([dict(res2[r]) for r in res2], index = res2.keys())
            res_df1.columns = ['clone_size_shannon_' + x for x in res_df1.columns]
            res_df2.columns = ['clone_degree_shannon_' + x for x in res_df2.columns]
            res_df = [res_df1, res_df2]
        return(res_df)

    def transfer_shannon(self, shannon_results, groupby, splitby = None):
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

    if return_dataframe:
        return(shannon_entropy(self, groupby, splitby, clone_key))
        logg.info(' finished', time=start)
    else:
        results = _shannon_entropy(self, groupby, splitby, clone_key)
        transfer_shannon(self, results, groupby, splitby)
        if self.__class__ == Dandelion:
            logg.info(' finished', time=start,
                deep=('updated `.metadata` with Shannon entropy\n'))
        if self.__class__ == AnnData:
            logg.info(' finished', time=start,
                deep=('updated `.obs` with Shannon entropy\n'))
