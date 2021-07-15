#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 14:01:32
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-07-15 23:01:41

import os
from collections import defaultdict, Iterable
import pandas as pd
import numpy as np
from subprocess import run
import re
from typing import Sequence, Tuple, Dict, Union
try:
    from typing import Literal
except ImportError:
    try:
        from typing_extensions import Literal
    except ImportError:

        class LiteralMeta(type):
            """LiteralMeta class."""

            def __getitem__(self, values):
                """Return Literal."""
                if not isinstance(values, tuple):
                    values = (values, )
                return type("Literal_", (Literal, ), dict(__args__=values))

        class Literal(metaclass=LiteralMeta):
            pass


class Tree(defaultdict):
    """Create a recursive defaultdict."""

    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value


def dict_from_table(meta: pd.DataFrame, columns: Tuple[str, str]) -> Dict:
    """
    Generate a dictionary from a dataframe.

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
    return (sample_dict)


def clean_nan_dict(d: Dict) -> Dict:
    """
    Remove nan from dictionary.

    Parameters
    ----------
    d : Dict
        dictionary

    Returns
    -------
        dictionary with no NAs.
    """
    return {k: v for k, v in d.items() if v is not np.nan}


def flatten(l: Sequence) -> Sequence:
    """
    Flatten a list-in-list-in-list.

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
    """
    cmd = ['makeblastdb', '-dbtype', 'nucl', '-parse_seqids', '-in', ref]
    run(cmd)


def bh(pvalues: np.array) -> np.array:
    """
    Compute the Benjamini-Hochberg FDR correction.

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
    return dataframe[key].dtype == str or dataframe[
        key].dtype == object or is_categorical(
            dataframe[key]) or dataframe[key].dtype == bool


def isGZIP(filename: str) -> bool:
    if filename.split('.')[-1] == 'gz':
        return True
    return False


def isBZIP(filename: str) -> bool:
    if filename.split('.')[-1] == 'pbz2':
        return True
    return False


def check_filepath(s,
                   filename_prefix: Union[None, str] = None,
                   endswith: Union[None, str] = None,
                   subdir: Union[None, str] = None):
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
                        out_ = out_ + os.path.basename(
                            file) + '/' + subdir + '/'
                    for x in os.listdir(out_):
                        if x.endswith(ends_with):
                            if str(x).split(
                                    ends_with)[0] == filename_pre + '_contig':
                                filePath = out_ + x
                            else:
                                continue
            else:
                continue
    return (filePath)


def check_fastapath(fasta, filename_prefix: Union[None, str] = None):
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
    return (filePath)


def change_file_location(data: Sequence,
                         filename_prefix: Union[None, Sequence, str] = None):
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


class FilterContigs:
    def __init__(self, data):
        self.data = data
        self.poor_qual = []
        self.h_doublet = []
        self.l_doublet = []
        self.drop_contig = []
        self.umi_adjustment = {}

    def run_scan(self, b, keep_highest_umi, umi_foldchange_cutoff, filter_poorqualitycontig, v_dict, d_dict, j_dict, c_dict):
        """Main workhorse of filter_contig."""
        h = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD']))]['sequence_id'])
        h_umi = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD']))]['duplicate_count']]
        if 'sequence_alignment' in self.data:
            h_seq = [x for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD']))]['sequence_alignment']]
        h_ccall = [x for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD']))]['c_call']]

        l = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]['sequence_id'])
        l_umi = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]['duplicate_count']]
        if 'sequence_alignment' in self.data:
            l_seq = [x for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]['sequence_alignment']]

        # split to productive vs non-productive
        h_p = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD'])) & (self.data['productive'].isin([True]))]['sequence_id'])
        h_umi_p = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD'])) & (self.data['productive'].isin([True]))]['duplicate_count']]
        h_np = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD'])) & (self.data['productive'].isin([False]))]['sequence_id'])
        h_umi_np = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD'])) & (self.data['productive'].isin([False]))]['duplicate_count']]
        l_p = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])) & (self.data['productive'].isin([True]))]['sequence_id'])
        l_umi_p = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])) & (self.data['productive'].isin([True]))]['duplicate_count']]
        l_np = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])) & (self.data['productive'].isin([False]))]['sequence_id'])
        l_umi_np = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])) & (self.data['productive'].isin([False]))]['duplicate_count']]

        # marking doublets defined by heavy chains
        if len(h) > 1:
            if 'sequence_alignment' in self.data:
                if len(list(set(h_seq))) == 1:
                    highest_umi_h = max(h_umi)
                    highest_umi_h_idx = [i for i, j in enumerate(h_umi) if j == highest_umi_h]
                    keep_index_h = highest_umi_h_idx[0]
                    self.drop_contig.append(h[:keep_index_h] + h[keep_index_h + 1:])
                    keep_hc_contig = h[keep_index_h]
                    self.data.at[keep_hc_contig, 'duplicate_count'] = int(np.sum(h_umi[:keep_index_h] + h_umi[keep_index_h + 1:]))
                    self.umi_adjustment.update({keep_hc_contig: int(np.sum(h_umi[:keep_index_h] + h_umi[keep_index_h + 1:]))})
                    h = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD']))]['sequence_id'])

                    # refresh
                    h_p = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD'])) & (self.data['productive'].isin([True]))]['sequence_id'])
                    h_umi_p = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD'])) & (self.data['productive'].isin([True]))]['duplicate_count']]
                    h_np = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD'])) & (self.data['productive'].isin([False]))]['sequence_id'])
                    h_umi_np = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD'])) & (self.data['productive'].isin([False]))]['duplicate_count']]
            if len(h_p) > 1:
                highest_umi_h = max(h_umi_p)
                highest_umi_idx = [i for i, j in enumerate(h_umi_p) if j == highest_umi_h]
                keep_index_h = highest_umi_idx[0]
                umi_test = [highest_umi_h / x < umi_foldchange_cutoff for x in h_umi_p[:keep_index_h] + h_umi_p[keep_index_h + 1:]]
                sum_umi = sum(h_umi_p)
                if 'IGHM' and 'IGHD' in h_ccall:
                    if all(cc_ == 'IGHM' or cc_ == 'IGHD' for cc_ in h_ccall):
                        pass
                else:
                    if len(highest_umi_idx) > 1:
                        self.h_doublet.append(b)
                    if sum_umi < 4:
                        self.h_doublet.append(b)
                    if any(umi_test):
                        self.h_doublet.append(b)
                    if len(highest_umi_idx) == 1:
                        other_umi_idx = [i for i, j in enumerate(h_umi_p)if j != highest_umi_h]
                        umi_test_ = [highest_umi_h / x >= umi_foldchange_cutoff for x in h_umi_p[:keep_index_h] + h_umi_p[keep_index_h + 1:]]
                        umi_test_dict = dict(zip(other_umi_idx, umi_test_))
                        for otherindex in umi_test_dict:
                            if umi_test_dict[otherindex]:
                                if keep_highest_umi:
                                    self.drop_contig.append(h_p[otherindex])
                                    # refresh
                                    h_p = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD'])) & (self.data['productive'].isin([True]))]['sequence_id'])

            # tries to keep most non-productive contigs but if there is a dominant contig, drop other contigs.
            if len(h_np) > 1:
                highest_umi_h = max(h_umi_np)
                highest_umi_idx = [i for i, j in enumerate(h_umi_np) if j == highest_umi_h]
                if len(highest_umi_idx) == 1:
                    keep_index_h = highest_umi_idx[0]
                    other_umi_idx = [i for i, j in enumerate(h_umi_np)if j != highest_umi_h]
                    umi_test_ = [highest_umi_h / x >= umi_foldchange_cutoff for x in h_umi_np[:keep_index_h] + h_umi_np[keep_index_h + 1:]]
                    umi_test_dict = dict(zip(other_umi_idx, umi_test_))
                    for otherindex in umi_test_dict:
                        if umi_test_dict[otherindex]:
                            self.drop_contig.append(h_np[otherindex])
                            # refresh
                            h_np = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD'])) & (self.data['productive'].isin([False]))]['sequence_id'])
        if len(l) > 1:
            if 'sequence_alignment' in self.data:
                if len(list(set(l_seq))) == 1:
                    highest_umi_l = max(l_umi)
                    highest_umi_l_idx = [i for i, j in enumerate(l_umi) if j == highest_umi_l]
                    keep_index_l = highest_umi_l_idx[0]
                    self.drop_contig.append(l[:keep_index_l] + l[keep_index_l + 1:])
                    keep_lc_contig = l[keep_index_l]
                    self.data.at[keep_lc_contig, 'duplicate_count'] = int(np.sum(l_umi[:keep_index_l] + l_umi[keep_index_l + 1:]))
                    self.umi_adjustment.update({keep_lc_contig: int(np.sum(l_umi[:keep_index_l] + l_umi[keep_index_l + 1:]))})
                    l = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]['sequence_id'])

                    # refresh
                    l_p = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])) & (self.data['productive'].isin([True]))]['sequence_id'])
                    l_umi_p = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])) & (self.data['productive'].isin([True]))]['duplicate_count']]
                    l_np = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])) & (self.data['productive'].isin([False]))]['sequence_id'])
                    l_umi_np = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])) & (self.data['productive'].isin([False]))]['duplicate_count']]
            if len(l_p) > 1:
                highest_umi_l = max(l_umi_p)
                highest_umi_l_idx = [i for i, j in enumerate(l_umi_p) if j == highest_umi_l]
                keep_index_l = highest_umi_l_idx[0]
                umi_test = [highest_umi_l / x < umi_foldchange_cutoff for x in h_umi_p[:keep_index_l] + h_umi_p[keep_index_l + 1:]]
                sum_umi = sum(l_umi_p)
                if len(highest_umi_l_idx) > 1:
                    self.l_doublet.append(b)
                if sum_umi < 4:
                    self.l_doublet.append(b)
                if any(umi_test):
                    self.l_doublet.append(b)
                if len(highest_umi_l_idx) == 1:
                    keep_index_l = highest_umi_l_idx[0]
                    other_umi_idx_l = [i for i, j in enumerate(l_umi_p) if j != highest_umi_l]
                    umi_test_l = [highest_umi_l / x >= umi_foldchange_cutoff for x in l_umi_p[:keep_index_l] + l_umi_p[keep_index_l + 1:]]
                    umi_test_dict_l = dict(zip(other_umi_idx_l, umi_test_l))
                    for otherindex in umi_test_dict_l:
                        if umi_test_dict_l[otherindex]:
                            if keep_highest_umi:
                                self.drop_contig.append(l_p[otherindex])
                                # refresh
                                l_p = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])) & (self.data['productive'].isin([True]))]['sequence_id'])

            # tries to keep most non-productive contigs but if there is a dominant contig, drop other contigs.
            if len(l_np) > 1:
                highest_umi_l = max(l_umi_np)
                highest_umi_l_idx = [i for i, j in enumerate(l_umi_np) if j == highest_umi_l]
                keep_index_l = highest_umi_l_idx[0]
                other_umi_idx_l = [i for i, j in enumerate(l_umi_np) if j != highest_umi_l]
                umi_test_l = [highest_umi_l / x >= umi_foldchange_cutoff for x in l_umi_np[:keep_index_l] + l_umi_np[keep_index_l + 1:]]
                if len(highest_umi_l_idx) == 1:
                    umi_test_dict_l = dict(zip(other_umi_idx_l, umi_test_l))
                    for otherindex in umi_test_dict_l:
                        if umi_test_dict_l[otherindex]:
                            if keep_highest_umi:
                                self.drop_contig.append(l_np[otherindex])
                                # refresh
                                l_np = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])) & (self.data['productive'].isin([False]))]['sequence_id'])

        # marking doublets defined by VJ chains
        if (len(h_p) == 1) & (len(l_p) > 1):
            self.l_doublet.append(b)

        # marking poor bcr quality, defined as those with only VJ chains, those
        # that were have conflicting assignment of locus and V(D)J v-, d-, j- and c- calls,
        # and also those that are missing j calls (to catch non-productive).
        if len(h) < 1:
            if filter_poorqualitycontig:
                self.poor_qual.append(b)
            self.drop_contig.append(l)
        if len(h_p) == 1:
            v = v_dict[h_p[0]]
            j = j_dict[h_p[0]]
            d = d_dict[h_p[0]]
            c = c_dict[h_p[0]]
            if pd.notnull(v) and v != '':
                if not re.search('IGH|TR[BD]', v):
                    if filter_poorqualitycontig:
                        self.poor_qual.append(b)
                    self.drop_contig.append(l_p)
                    self.drop_contig.append(h_p)
            if pd.notnull(d) and d != '':
                if not re.search('IGH|TR[BD]', d):
                    if filter_poorqualitycontig:
                        self.poor_qual.append(b)
                    self.drop_contig.append(l_p)
                    self.drop_contig.append(h_p)
            if pd.notnull(j) and j != '':
                if not re.search('IGH|TR[BD]', j):
                    if filter_poorqualitycontig:
                        self.poor_qual.append(b)
                    self.drop_contig.append(l_p)
                    self.drop_contig.append(h_p)
            if pd.notnull(c) and c != '':
                if not re.search('IGH|TR[BD]', c):
                    if filter_poorqualitycontig:
                        self.poor_qual.append(b)
                    self.drop_contig.append(l_p)
                    self.drop_contig.append(h_p)

            if pd.notnull(j) and j != '':
                if pd.notnull(v) and v != '':
                    if (re.search('IGH', v) and not re.search('IGH', j)) or (re.search('IGH', j) and not re.search('IGH', v)):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(l_p)
                        self.drop_contig.append(h_p)
                    elif (re.search('TRB', v) and not re.search('TRB', j)) or (re.search('TRB', j) and not re.search('TRB', v)):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(l_p)
                        self.drop_contig.append(h_p)
                    elif (re.search('TRD', v) and not re.search('TRD', j)) or (re.search('TRD', j) and not re.search('TRD', v)):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(l_p)
                        self.drop_contig.append(h_p)

                if pd.notnull(d) and d != '':
                    if (re.search('IGH', d) and not re.search('IGH', j)) or (re.search('IGH', d) and not re.search('IGH', v)):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(l_p)
                        self.drop_contig.append(h_p)
                    elif (re.search('TRB', d) and not re.search('TRB', j)) or (re.search('TRB', d) and not re.search('TRB', v)):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(l_p)
                        self.drop_contig.append(h_p)
                    elif (re.search('TRD', d) and not re.search('TRD', j)) or (re.search('TRD', d) and not re.search('TRD', v)):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(l_p)
                        self.drop_contig.append(h_p)
            else:
                if filter_poorqualitycontig:
                    self.poor_qual.append(b)
                self.drop_contig.append(l_p)
                self.drop_contig.append(h_p)

        if len(h_p) > 1:
            for hx in h_p:
                v = v_dict[hx]
                d = d_dict[hx]
                j = j_dict[hx]
                c = c_dict[hx]
                if pd.notnull(v) and v != '':
                    if not re.search('IGH|TR[BD]', v):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(hx)
                if pd.notnull(d) and d != '':
                    if not re.search('IGH|TR[BD]', d):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(hx)
                if pd.notnull(j) and j != '':
                    if not re.search('IGH|TR[BD]', j):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(hx)
                if pd.notnull(c) and c != '':
                    if not re.search('IGH|TR[BD]', c):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(hx)
                if pd.notnull(j) and j != '':
                    if pd.notnull(v) and v != '':
                        if (re.search('IGH', v) and not re.search('IGH', j)) or (re.search('IGH', j) and not re.search('IGH', v)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(hx)
                        elif (re.search('TRB', v) and not re.search('TRB', j)) or (re.search('TRB', j) and not re.search('TRB', v)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(hx)
                        elif (re.search('TRD', v) and not re.search('TRD', j)) or (re.search('TRD', j) and not re.search('TRD', v)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(hx)
                    if pd.notnull(d) and d != '':
                        if (re.search('IGH', d) and not re.search('IGH', j)) or (re.search('IGH', d) and not re.search('IGH', v)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(hx)
                        elif (re.search('TRB', d) and not re.search('TRB', j)) or (re.search('TRB', d) and not re.search('TRB', v)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(hx)
                        elif (re.search('TRD', d) and not re.search('TRD', j)) or (re.search('TRD', d) and not re.search('TRD', v)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(hx)
                else:
                    if filter_poorqualitycontig:
                        self.poor_qual.append(b)
                    self.drop_contig.append(hx)

        if len(h_np) > 0:
            for hx in h_np:
                v = v_dict[hx]
                d = d_dict[hx]
                j = j_dict[hx]
                c = c_dict[hx]
                if pd.notnull(v) and v != '':
                    if not re.search('IGH|TR[BD]', v):
                        self.drop_contig.append(hx)
                if pd.notnull(d) and d != '':
                    if not re.search('IGH|TR[BD]', d):
                        self.drop_contig.append(hx)
                if pd.notnull(j) and j != '':
                    if not re.search('IGH|TR[BD]', j):
                        self.drop_contig.append(hx)
                if pd.notnull(c) and c != '':
                    if not re.search('IGH|TR[BD]', c):
                        self.drop_contig.append(hx)

                if pd.notnull(j) and j != '':
                    if pd.notnull(v) and v != '':
                        if (re.search('IGH', v) and not re.search('IGH', j)) or (re.search('IGH', j) and not re.search('IGH', v)):
                            self.drop_contig.append(hx)
                        elif (re.search('TRB', v) and not re.search('TRB', j)) or (re.search('TRB', j) and not re.search('TRB', v)):
                            self.drop_contig.append(hx)
                        elif (re.search('TRD', v) and not re.search('TRD', j)) or (re.search('TRD', j) and not re.search('TRD', v)):
                            self.drop_contig.append(hx)
                    if pd.notnull(d) and d != '':
                        if (re.search('IGH', d) and not re.search('IGH', j)) or (re.search('IGH', d) and not re.search('IGH', v)):
                            self.drop_contig.append(hx)
                        elif (re.search('TRB', d) and not re.search('TRB', j)) or (re.search('TRB', d) and not re.search('TRB', v)):
                            self.drop_contig.append(hx)
                        elif (re.search('TRD', d) and not re.search('TRD', j)) or (re.search('TRD', d) and not re.search('TRD', v)):
                            self.drop_contig.append(hx)
                else:
                    self.drop_contig.append(hx)
        if len(l_p) > 0:
            for lx in l_p:
                v = v_dict[lx]
                j = j_dict[lx]
                c = c_dict[lx]
                if pd.notnull(v) and v != '':
                    if re.search('IGH|TR[BD]', v):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(lx)
                if pd.notnull(j) and j != '':
                    if re.search('IGH|TR[BD]', j):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(lx)
                if pd.notnull(c) and c != '':
                    if re.search('IGH|TR[BD]', c):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(lx)

                if pd.notnull(j) and j != '':
                    if pd.notnull(v) and v != '':
                        if (re.search('IGK', v) and not re.search('IGK', j)) or (re.search('IGK', j) and not re.search('IGK', v)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                        if (re.search('IGL', v) and not re.search('IGL', j)) or (re.search('IGL', j) and not re.search('IGL', v)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                        elif (re.search('TRA', v) and not re.search('TRA', j)) or (re.search('TRA', j) and not re.search('TRA', v)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                        elif (re.search('TRG', v) and not re.search('TRG', j)) or (re.search('TRG', j) and not re.search('TRG', v)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                else:
                    if filter_poorqualitycontig:
                        self.poor_qual.append(b)
                    self.drop_contig.append(lx)

        if len(l_np) > 0:
            for lx in l_np:
                v = v_dict[lx]
                j = j_dict[lx]
                c = c_dict[lx]
                if pd.notnull(v) and v != '':
                    if re.search('IGH|TR[BD]', v):
                        self.drop_contig.append(lx)
                if pd.notnull(j) and j != '':
                    if re.search('IGH|TR[BD]', j):
                        self.drop_contig.append(lx)
                if pd.notnull(c) and c != '':
                    if re.search('IGH|TR[BD]', c):
                        self.drop_contig.append(lx)

                if pd.notnull(j) and j != '':
                    if pd.notnull(v) and v != '':
                        if (re.search('IGK', v) and not re.search('IGK', j)) or (re.search('IGK', j) and not re.search('IGK', v)):
                            self.drop_contig.append(lx)
                        if (re.search('IGL', v) and not re.search('IGL', j)) or (re.search('IGL', j) and not re.search('IGL', v)):
                            self.drop_contig.append(lx)
                        elif (re.search('TRA', v) and not re.search('TRA', j)) or (re.search('TRA', j) and not re.search('TRA', v)):
                            self.drop_contig.append(lx)
                        elif (re.search('TRG', v) and not re.search('TRG', j)) or (re.search('TRG', j) and not re.search('TRG', v)):
                            self.drop_contig.append(lx)
                else:
                    self.drop_contig.append(lx)

    def run_scan_lite(self, b, v_dict, d_dict, j_dict, c_dict):
        """A 'lite' version of filter_contig where it does not impose the 1 VDJ + 1 VJ filtering."""
        h = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD']))]['sequence_id'])
        h_umi = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD']))]['duplicate_count']]
        if 'sequence_alignment' in self.data:
            h_seq = [x for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD']))]['sequence_alignment']]

        l = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]['sequence_id'])
        l_umi = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]['duplicate_count']]
        if 'sequence_alignment' in self.data:
            l_seq = [x for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]['sequence_alignment']]

        if len(h) > 1:
            if 'sequence_alignment' in self.data:
                if len(list(set(h_seq))) == 1:
                    highest_umi_h = max(h_umi)
                    highest_umi_h_idx = [i for i, j in enumerate(h_umi) if j == highest_umi_h]
                    keep_index_h = highest_umi_h_idx[0]
                    self.drop_contig.append(h[:keep_index_h] + h[keep_index_h + 1:])
                    keep_hc_contig = h[keep_index_h]
                    self.data.at[keep_hc_contig, 'duplicate_count'] = int(np.sum(h_umi[:keep_index_h] + h_umi[keep_index_h + 1:]))
                    self.umi_adjustment.update({keep_hc_contig: int(np.sum(h_umi[:keep_index_h] + h_umi[keep_index_h + 1:]))})
                    h = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD']))]['sequence_id'])
                    h_umi = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGH', 'TRB', 'TRD']))]['duplicate_count']]

        if len(l) > 1:
            if 'sequence_alignment' in self.data:
                if len(list(set(l_seq))) == 1:
                    highest_umi_l = max(l_umi)
                    highest_umi_l_idx = [i for i, j in enumerate(l_umi) if j == highest_umi_l]
                    keep_index_l = highest_umi_l_idx[0]
                    self.drop_contig.append(l[:keep_index_l] + l[keep_index_l + 1:])
                    keep_hc_contig = l[keep_index_l]
                    self.data.at[keep_hc_contig, 'duplicate_count'] = int(np.sum(l_umi[:keep_index_l] + l_umi[keep_index_l + 1:]))
                    self.umi_adjustment.update({keep_lc_contig: int(np.sum(l_umi[:keep_index_l] + l_umi[keep_index_l + 1:]))})
                    l = list(self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]['sequence_id'])
                    l_umi = [int(x) for x in self.data[(self.data['cell_id'].isin([b])) & (self.data['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]['duplicate_count']]

        if len(h) > 0:
            for hx in h:
                v = v_dict[hx]
                d = d_dict[hx]
                j = j_dict[hx]
                c = c_dict[hx]
                if pd.notnull(v) and v != '':
                    if not re.search('IGH|TR[BD]', v):
                        self.drop_contig.append(hx)
                if pd.notnull(d) and d != '':
                    if not re.search('IGH|TR[BD]', d):
                        self.drop_contig.append(hx)
                if pd.notnull(j) and j != '':
                    if not re.search('IGH|TR[BD]', j):
                        self.drop_contig.append(hx)
                if pd.notnull(c) and c != '':
                    if not re.search('IGH|TR[BD]', c):
                        self.drop_contig.append(hx)

                if pd.notnull(j) and j != '':
                    if pd.notnull(v) and v != '':
                        if (re.search('IGH', v) and not re.search('IGH', j)) or (re.search('IGH', j) and not re.search('IGH', v)):
                            self.drop_contig.append(hx)
                        elif (re.search('TRB', v) and not re.search('TRB', j)) or (re.search('TRB', j) and not re.search('TRB', v)):
                            self.drop_contig.append(hx)
                        elif (re.search('TRD', v) and not re.search('TRD', j)) or (re.search('TRD', j) and not re.search('TRD', v)):
                            self.drop_contig.append(hx)
                    if pd.notnull(d) and d != '':
                        if (re.search('IGH', d) and not re.search('IGH', j)) or (re.search('IGH', d) and not re.search('IGH', v)):
                            self.drop_contig.append(hx)
                        elif (re.search('TRB', d) and not re.search('TRB', j)) or (re.search('TRB', d) and not re.search('TRB', v)):
                            self.drop_contig.append(hx)
                        elif (re.search('TRD', d) and not re.search('TRD', j)) or (re.search('TRD', d) and not re.search('TRD', v)):
                            self.drop_contig.append(hx)
                else:
                    self.drop_contig.append(hx)
        if len(l) > 0:
            for lx in l:
                v = v_dict[lx]
                j = j_dict[lx]
                c = c_dict[lx]
                if pd.notnull(v) and v != '':
                    if re.search('IGH|TR[BD]', v):
                        self.drop_contig.append(lx)
                if pd.notnull(j) and j != '':
                    if re.search('IGH|TR[BD]', j):
                        self.drop_contig.append(lx)
                if pd.notnull(c) and c != '':
                    if re.search('IGH|TR[BD]', c):
                        self.drop_contig.append(lx)

                if pd.notnull(j) and j != '':
                    if pd.notnull(v) and v != '':
                        if (re.search('IGK', v) and not re.search('IGK', j)) or (re.search('IGK', j) and not re.search('IGK', v)):
                            self.drop_contig.append(lx)
                        if (re.search('IGL', v) and not re.search('IGL', j)) or (re.search('IGL', j) and not re.search('IGL', v)):
                            self.drop_contig.append(lx)
                        elif (re.search('TRA', v) and not re.search('TRA', j)) or (re.search('TRA', j) and not re.search('TRA', v)):
                            self.drop_contig.append(lx)
                        elif (re.search('TRG', v) and not re.search('TRG', j)) or (re.search('TRG', j) and not re.search('TRG', v)):
                            self.drop_contig.append(lx)
                else:
                    self.drop_contig.append(lx)


def cmp(a, b):
    """Python2.x cmp function."""
    return (a > b) - (a < b)


def cmp_str_emptylast(s1, s2):
    """Help sort empty string to last."""
    if not s1 or not s2:
        return bool(s2) - bool(s1)

    return cmp(s1, s2)


def cmp_to_key(mycmp):
    """Convert a cmp= function into a key= function."""
    class K:
        def __init__(self, obj, *args):
            self.obj = obj

        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0

        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0

        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0

        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0

        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0

        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0

    return K
