#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 14:01:32
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-06-15 10:05:47

import os
from collections import defaultdict, Iterable
import pandas as pd
import numpy as np
from subprocess import run
from typing import Sequence, Tuple, Dict, Union
try:
    from typing import Literal
except ImportError:
    try:
        from typing_extensions import Literal
    except ImportError:
        class LiteralMeta(type):
            def __getitem__(cls, values):
                if not isinstance(values, tuple):
                    values = (values,)
                return type("Literal_", (Literal,), dict(__args__=values))

        class Literal(metaclass=LiteralMeta):
            pass


class Tree(defaultdict):
    '''
    Create a recursive defaultdict
    '''

    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value


def dict_from_table(meta: pd.DataFrame, columns: Tuple[str, str]) -> Dict:
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


def clean_nan_dict(d: Dict) -> Dict:
    """
    Parameters
    ----------
    d : Dict
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


def flatten(l: Sequence) -> Sequence:
    """
    Parameters
    ----------
    l : Sequence
        a list-in-list list

    Returns
    -------
        a flattened list.
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


def makeblastdb(ref: str):
    """
    Run makeblastdb on constant region fasta file.

    Wrapper for makeblastdb.

    Parameters
    ----------
    ref : str
        constant region fasta file
    Returns
    -------


    """

    cmd = ['makeblastdb',
           '-dbtype', 'nucl',
           '-parse_seqids',
           '-in', ref]
    run(cmd)


def bh(pvalues: np.array) -> np.array:
    """
    Computes the Benjamini-Hochberg FDR correction.

    Parameters
    ----------
        pvalues : np.array
            array of p-values to correct
    Returns
    -------
        np.array of corrected p-values
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
        new_values.append((n / rank) * pvalue)
    for i in range(0, int(n) - 1):
        if new_values[i] < new_values[i + 1]:
            new_values[i + 1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues


def is_categorical(array_like) -> bool:
    return array_like.dtype.name == 'category'


def type_check(dataframe, key) -> bool:
    return dataframe[key].dtype == str or dataframe[key].dtype == object or is_categorical(dataframe[key]) or dataframe[key].dtype == bool


def isGZIP(filename: str) -> bool:
    if filename.split('.')[-1] == 'gz':
        return True
    return False


def isBZIP(filename: str) -> bool:
    if filename.split('.')[-1] == 'pbz2':
        return True
    return False


def check_filepath(s, filename_prefix = None, endswith = None, subdir = None):
    if filename_prefix is None:
        filename_pre = 'filtered'
    else:
        filename_pre = filename_prefix

    if endswith is None:
        ends_with = ''
    else:
        ends_with = endswith

    filePath = None
    if os.path.isfile(str(s)) and str(s).endswith(ends_with):
        filePath = s
    elif os.path.isdir(str(s)):
        files = os.listdir(s)
        for file in files:
            out_ = s.rstrip('/') + '/'
            if os.path.isdir(out_ + os.path.basename(file)):
                if file == 'dandelion':
                    if subdir is None:
                        out_ = out_ + os.path.basename(file) + '/'
                    else:
                        out_ = out_ + os.path.basename(file) + '/' + subdir + '/'
                    for x in os.listdir(out_):
                        if x.endswith(ends_with):
                            if str(x).split(ends_with)[0] == filename_pre + '_contig':
                                filePath = out_ + x
                            else:
                                continue
            else:
                continue
    return(filePath)


def check_fastapath(fasta, filename_prefix = None):
    if filename_prefix is None:
        filename_pre = 'filtered'
    else:
        filename_pre = filename_prefix

    filePath = None
    if os.path.isfile(str(fasta)) and str(fasta).endswith('.fasta'):
        filePath = fasta
    elif os.path.isdir(str(fasta)):
        files = os.listdir(fasta)
        for file in files:
            out_ = fasta.rstrip('/') + '/'
            if str(file).endswith(".fasta"):
                if str(file).split('.fasta')[0] == filename_pre + '_contig':
                    filePath = out_ + os.path.basename(file)
                else:
                    continue
    return(filePath)


def change_file_location(data: Sequence, filename_prefix: Union[None, Sequence, str] = None):
    """
    Move file from tmp folder to dandelion folder.

    Only used for TCR data.

    Parameters
    ----------
    data : Sequence
        list of data folders containing the .tsv files. if provided as a single string, it will first be converted to a
        list; this allows for the function to be run on single/multiple samples.
    filename_prefix : str, optional
        list of prefixes of file names preceding '_contig'. None defaults to 'filtered'.

    Returns
    -------
    Individual V(D)J data files with v_call_genotyped column containing reassigned heavy chain v calls
    """
    fileformat = 'blast'
    if type(data) is not list:
        data = [data]
    if type(filename_prefix) is not list:
        filename_prefix = [filename_prefix]
    if all(t is None for t in filename_prefix):
        filename_prefix = [None for d in data]

    informat_dict = {
        'changeo': '_igblast_db-pass.tsv',
        'blast': '_igblast_db-pass.tsv',
        'airr': '_igblast_gap.tsv'
    }

    filePath = None

    for i in range(0, len(data)):
        filePath = check_filepath(data[i],
                                  filename_prefix=filename_prefix[i],
                                  endswith=informat_dict[fileformat],
                                  subdir='tmp')
        if filePath is None:
            raise OSError(
                'Path to .tsv file for {} is unknown. '.format(data[i]) +
                'Please specify path to reannotated .tsv file or folder containing reannotated .tsv file.'
            )

        cmd = ['rsync', '-azvh', filePath, filePath.rsplit('/', 2)[0]]

        run(cmd)
