#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 14:01:32
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-02-11 12:24:47

import os
from collections import defaultdict, Iterable
import pandas as pd
import numpy as np
from subprocess import run


class Tree(defaultdict):
    '''
    Create a recursive defaultdict
    '''

    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value


def dict_from_table(meta, columns):
    """
    Generates a dictionary from a dataframe
    Parameters
    ----------
    meta
        pandas dataframe or file path
    columns
        column names in dataframe

    Returns
    -------
    dictionary
    """
    if (isinstance(meta, pd.DataFrame)) & (columns is not None):
        meta_ = meta
        if len(columns) == 2:
            sample_dict = dict(zip(meta_[columns[0]], meta_[columns[1]]))
    elif (os.path.isfile(str(meta))) & (columns is not None):
        meta_ = pd.read_csv(meta, sep='\t', dtype='object')
        if len(columns) == 2:
            sample_dict = dict(zip(meta_[columns[0]], meta_[columns[1]]))

    sample_dict = clean_nan_dict(sample_dict)
    return(sample_dict)


def clean_nan_dict(d):
    """
    Parameters
    ----------
    d
        dictionary

    Returns
    -------
    dictionary with no NAs.
    """

    return {
        k: v
        for k, v in d.items()
        if v is not np.nan
    }


def flatten(l):
    """
    Parameters
    ----------
    l
        list

    Returns
    -------
    a flattened list.
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


def makeblastdb(ref):
    """
    Runs makeblastdb.

    Parameters
    ----------
    ref
        fasta file

    """

    cmd = ['makeblastdb',
           '-dbtype', 'nucl',
           '-parse_seqids',
           '-in', ref]
    run(cmd)


def bh(pvalues):
    """
    Computes the Benjamini-Hochberg FDR correction.

    Input:
        * pvals - vector of p-values to correct
    """
    n = int(pvalues.shape[0])
    new_pvalues = np.empty(n)
    values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues


def is_categorical(array_like):
    return array_like.dtype.name == 'category'


def type_check(dataframe, key):
    return dataframe[key].dtype == str or dataframe[key].dtype == object or is_categorical(dataframe[key]) or dataframe[key].dtype == bool


def isGZIP(filename):
    if filename.split('.')[-1] == 'gz':
        return True
    return False


def isBZIP(filename):
    if filename.split('.')[-1] == 'pbz2':
        return True
    return False
