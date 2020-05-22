#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 14:01:32
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-05-22 14:58:08

import sys
import os
from collections import defaultdict, Iterable
import pandas as pd
import numpy as np
from subprocess import run
from tqdm import tqdm
import re

class Tree(defaultdict):
    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value

def fasta_iterator(fh):
    while True:
        line = fh.readline()
        if line.startswith('>'):
            break
    while True:
        header = line[1:-1].rstrip()
        sequence = fh.readline().rstrip()
        while True:
            line = fh.readline()
            if not line:
                break
            if line.startswith('>'):
                break
            sequence += line.rstrip()
        yield(header, sequence)
        if not line:
            return

def Write_output(out, file):
    fh = open(file, "a")
    fh.write(out)
    fh.close()
    return()

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
    sample_dict
        dictionary
    """
    if (isinstance(meta, pd.DataFrame)) & (columns is not None):
        meta_ = meta
        if len(columns) == 2:
            sample_dict = dict(zip(meta_[columns[0]], meta_[columns[1]]))
    elif (os.path.isfile(str(meta))) & (columns is not None):
        meta_ = pd.read_csv(meta, sep = '\t', dtype = 'object')
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
        k:v
        for k, v in d.items()
        if v is not np.nan
    }

def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el

def makeblastdb(ref):
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
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
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

def extract(d, keys):
    return(dict((k, d[k]) for k in keys if k in d))

def load_data(obj):
    if os.path.isfile(str(obj)):
        try:
            obj_ = pd.read_csv(obj, sep = '\t', dtype = 'object')
        except FileNotFoundError as e:
            print(e)
    elif isinstance(obj, pd.DataFrame):
            obj_ = obj.copy()
    else:
        raise TypeError("Input is not of <class 'pandas.core.frame.DataFrame'>")
        
    if 'sequence_id' in obj_.columns:
        obj_.set_index('sequence_id', drop = False, inplace = True)
    else:
        raise KeyError("'sequence_id' not found in columns of input")
        
    return(obj_)

def convert_preprocessed_tcr_10x(file, prefix = None, save = None):
        
    cr_annot = pd.read_csv(file)
    if prefix is not None:
        cr_annot['index']=[prefix+'_'+i for i in cr_annot['contig_id']]
    else:
        cr_annot['index']=[i for i in cr_annot['contig_id']]
    cr_annot.set_index('index', inplace = True)
    
    ddl_annot = pd.read_csv("{}/dandelion/data/tmp/{}_igblast.tsv".format(os.path.dirname(file), os.path.basename(file).split('_annotations.csv')[0]), sep = '\t')
    ddl_annot.set_index('sequence_id', inplace = True, drop = False)
    
    for i in tqdm(ddl_annot.index, desc = 'Processing data '):
        v = ddl_annot.loc[i, 'v_call']
        d = ddl_annot.loc[i, 'd_call']
        j = ddl_annot.loc[i, 'j_call']
        if v is not np.nan:
            ddl_annot.loc[i, 'v_call_igblast'] = ','.join(list(set(re.sub('[*][0-9][0-9]', '', v).split(','))))
            v_ = list(set(re.sub('[*][0-9][0-9]', '', v).split(',')))
            if re.match('IGH', ','.join(v_)):
                ddl_annot.loc[i, 'locus'] = 'IGH'
            if re.match('IGK', ','.join(v_)):
                ddl_annot.loc[i, 'locus'] = 'IGK'
            if re.match('IGL', ','.join(v_)):
                ddl_annot.loc[i, 'locus'] = 'IGL'
            if len(v_) > 1:
                ddl_annot.loc[i, 'locus'] = 'Multi'
        if d is not np.nan:
            ddl_annot.loc[i, 'd_call_igblast'] = ','.join(list(set(re.sub('[*][0-9][0-9]', '', d).split(','))))
            d_ = list(set(re.sub('[*][0-9][0-9]', '', d).split(',')))
            if re.match('IGH', ','.join(d_)):
                ddl_annot.loc[i, 'locus'] = 'IGH'
            if re.match('IGK', ','.join(d_)):
                ddl_annot.loc[i, 'locus'] = 'IGK'
            if re.match('IGL', ','.join(d_)):
                ddl_annot.loc[i, 'locus'] = 'IGL'
            if len(d_) > 1:
                ddl_annot.loc[i, 'locus'] = 'Multi'
        if j is not np.nan:
            ddl_annot.loc[i, 'j_call_igblast'] = ','.join(list(set(re.sub('[*][0-9][0-9]', '', j).split(','))))
            j_ = list(set(re.sub('[*][0-9][0-9]', '', j).split(',')))
            if re.match('IGH', ','.join(j_)):
                ddl_annot.loc[i, 'locus'] = 'IGH'
            if re.match('IGK', ','.join(j_)):
                ddl_annot.loc[i, 'locus'] = 'IGK'
            if re.match('IGL', ','.join(j_)):
                ddl_annot.loc[i, 'locus'] = 'IGL'
            if len(j_) > 1:
                ddl_annot.loc[i, 'locus'] = 'Multi'
    ddl_annot['cell_id'] = [c.split('_contig')[0].split('-')[0] for c in ddl_annot['sequence_id']]

    cellrangermap = {
        'cell_id':'barcode',
        'sequence_id':'contig_id',
        'locus':'chain',
        'v_call_igblast':'v_gene',
        'd_call_igblast':'d_gene',
        'j_call_igblast':'j_gene',
        'productive':'productive',
        'junction_aa':'cdr3',
        'junction':'cdr3_nt'}

    for i in tqdm(cr_annot.index, desc = 'Updating data'):
        for key, value in cellrangermap.items():        
            if cr_annot.loc[i, 'chain'] not in ['IGH', 'IGK', 'IGL', None]:
                cr_annot.loc[i, value] = ddl_annot.loc[i, key]
    if save is not None:
        cr_annot.to_csv(save, index = False)
    else:
        cr_annot.to_csv("{}/dandelion/data/{}".format(os.path.dirname(file), os.path.basename(file)), index = False)