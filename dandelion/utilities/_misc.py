#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 14:01:32
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-06-01 22:20:18

import sys
import os
from collections import defaultdict, Iterable
import pandas as pd
import numpy as np
from subprocess import run
from tqdm import tqdm
import re
import gzip
import pickle as pickle
import copy

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
        raise TypeError("Either input is not of <class 'pandas.core.frame.DataFrame'> or file does not exist.")

    if 'sequence_id' in obj_.columns:
        obj_.set_index('sequence_id', drop = False, inplace = True)
    else:
        raise KeyError("'sequence_id' not found in columns of input")

    return(obj_)

def initialize_metadata(self, retrieve = None, clones_sep = None):
    # a quick way to retrieve the 'meta data' for each cell that transfer into obs lot in scanpy later
    dat = load_data(self.data)
    for x in ['cell_id', 'locus', 'clone_id']:
        if x not in dat.columns:
            raise KeyError ("Please check your object. %s is not in the columns of input data." % x)
    dat_h = dat[dat['locus'] == 'IGH']
    dat_l = dat[~(dat['locus'] == 'IGH')]
    clone_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h['clone_id'])))
    clone_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l['clone_id'])))
    metadata = pd.DataFrame.from_dict(clone_h, orient = 'index', columns = ['cell_id', 'heavy'])
    metadata.set_index('cell_id', inplace = True)
    light_clone_tree = Tree()
    for key, value in clone_l.items():
        k, v = value
        light_clone_tree[k][key] = v
    light_clone_tree2 = Tree()
    for g in light_clone_tree:
        second_key = []
        for k2 in light_clone_tree[g].keys():
            second_key.append(k2)
        second_key = list(set(second_key))
        second_key_dict = dict(zip(second_key, range(0,len(second_key))))
        for key, value in light_clone_tree[g].items():
            light_clone_tree2[g][second_key_dict[key]] = value
    metadata['light'] = pd.Series(light_clone_tree2)
    tmp_dat = metadata['light'].apply(pd.Series)
    tmp_dat.columns = ['light_' + str(c) for c in tmp_dat.columns]
    metadata = metadata.merge(tmp_dat, left_index = True, right_index = True)
    metadata = metadata[['heavy'] + [str(c) for c in tmp_dat.columns]]
    clones_list = {}
    for x in metadata.index:
        cl = list(set(list(metadata.loc[x, :])))
        cl = sorted([y for y in cl if str(y) != 'nan'])
        if len(cl) > 1:
            cl = cl[1:]
        clones_list[x] = ','.join(cl)
    metadata['clone_id'] = pd.Series(clones_list)
    metadata = metadata[['clone_id']]
    if clones_sep is None:
        scb = (0, '_')
    else:
        scb = (clones_sep[0], clones_sep[1])
    group = []
    for x in metadata['clone_id']:
        if scb[1] not in x:
            warnings.warn(UserWarning("\n\nClones do not contain '{}' as separator. Will not split the clone.\n".format(scb[1])))
            group.append(x)
        else:
            group.append(x.split(scb[1])[scb[0]])
    metadata['clone_group_id'] = group
    # 4) Heavy chain V, J and Isotype (C)
    iso_dict = dict(zip(dat_h['cell_id'], dat_h['c_call']))
    if 'v_call_genotyped' in dat_h.columns:
        v_dict = dict(zip(dat_h['cell_id'], dat_h['v_call_genotyped']))
    else:
        v_dict = dict(zip(dat_h['cell_id'], dat_h['v_call']))
    j_dict = dict(zip(dat_h['cell_id'], dat_h['j_call']))
    metadata['isotype'] = pd.Series(iso_dict)
    metadata['heavychain_v'] = pd.Series(v_dict)
    metadata['heavychain_j'] = pd.Series(j_dict)
    
    # 5) light chain V, J and C
    lc_dict = dict(zip(dat_l['cell_id'], dat_l['c_call']))
    if 'v_call_genotyped' in dat_l.columns:
        vl_dict = dict(zip(dat_l['cell_id'], dat_l['v_call_genotyped']))
    else:
        vl_dict = dict(zip(dat_l['cell_id'], dat_l['v_call']))
    jl_dict = dict(zip(dat_l['cell_id'], dat_l['j_call']))
    metadata['lightchain'] = pd.Series(lc_dict)
    metadata['lightchain_v'] = pd.Series(vl_dict)
    metadata['lightchain_j'] = pd.Series(jl_dict)

    hv =[re.sub('[*][0-9][0-9]', '', str(v)) for v in metadata['heavychain_v']]
    hv = [','.join(list(set(v.split(',')))) for v in hv]
    hj =[re.sub('[*][0-9][0-9]', '', str(j)) for j in metadata['heavychain_j']]
    hj = [','.join(list(set(j.split(',')))) for j in hj]
    lv =[re.sub('[*][0-9][0-9]', '', str(v)) for v in metadata['lightchain_v']]
    lv = [','.join(list(set(v.split(',')))) for v in lv]
    lj =[re.sub('[*][0-9][0-9]', '', str(j)) for j in metadata['lightchain_j']]
    lj = [','.join(list(set(j.split(',')))) for j in lj]

    metadata['heavychain_v_gene'] = hv
    metadata['heavychain_j_gene'] = hj
    metadata['lightchain_v_gene'] = lv
    metadata['lightchain_j_gene'] = lj

    metadata['heavy_multi'] = False
    metadata['light_multi'] = False
    for i in metadata.index:
        hv_ = metadata.loc[i, 'heavychain_v_gene'].split(',')
        hj_ = metadata.loc[i, 'heavychain_j_gene'].split(',')
        lv_ = metadata.loc[i, 'lightchain_v_gene'].split(',')
        lj_ = metadata.loc[i, 'lightchain_j_gene'].split(',')
        if any(x > 1 for x in [len(hv_), len(hj_)]):
            metadata.loc[i, 'heavy_multi'] = True
        if any(x > 1 for x in [len(lv_), len(lj_)]):
            metadata.loc[i, 'light_multi'] = True

    # 2) whether or not chains are productive
    productive_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h['productive'])))
    productive_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l['productive'])))
    metadata_pro = pd.DataFrame.from_dict(productive_h, orient = 'index', columns = ['cell_id', 'heavy_pro'])
    metadata_pro.set_index('cell_id', inplace = True)
    light_productive_tree = Tree()
    for key, value in productive_l.items():
        k, v = value
        light_productive_tree[k][key] = v
    light_productive_tree2 = Tree()
    for g in light_productive_tree:
        second_key = []
        for k2 in light_productive_tree[g].keys():
            second_key.append(k2)
        second_key = list(set(second_key))
        second_key_dict = dict(zip(second_key, range(0,len(second_key))))
        for key, value in light_productive_tree[g].items():
            light_productive_tree2[g][second_key_dict[key]] = value
    metadata_pro['light_pro'] = pd.Series(light_productive_tree2)
    tmp_dat = metadata_pro['light_pro'].apply(pd.Series)
    tmp_dat.columns = ['light_pro_' + str(c) for c in tmp_dat.columns]
    metadata_pro = metadata_pro.merge(tmp_dat, left_index = True, right_index = True)
    metadata_pro = metadata_pro[['heavy_pro'] + [str(c) for c in tmp_dat.columns]]
    productive_list = {}
    for x in metadata_pro.index:
        cl = list(set(list(metadata_pro.loc[x, :])))
        cl = sorted([y for y in cl if str(y) != 'nan'])
        if len(cl) > 1:
            cl = cl[1:]
        productive_list[x] = ','.join(cl)
    metadata['productive'] = pd.Series(productive_list)

    # return this in this order
    metadata = metadata[['clone_id', 'clone_group_id', 'isotype', 'lightchain', 'productive', 'heavy_multi', 'light_multi', 'heavychain_v_gene', 'lightchain_v_gene', 'heavychain_j_gene', 'lightchain_j_gene', 'heavychain_v', 'lightchain_v', 'heavychain_j', 'lightchain_j']]
    
    # new function to retrieve non-standard columns
    if retrieve is not None:
        if retrieve in dat.columns:
            # non-standard options to transfer to metadata
            h_retrieve_dict = dict(zip(dat_h['cell_id'], dat_h[retrieve]))
            l_retrieve_dict = dict(zip(dat_l['cell_id'], dat_l[retrieve]))
            metadata[str(retrieve) + '_heavychain'] = pd.Series(h_retrieve_dict)
            metadata[str(retrieve) + '_lightchain'] = pd.Series(l_retrieve_dict)
        else:
            raise KeyError('Unknown column : \'%s\' to retrieve.' % retrieve)

    if self.metadata is None:
        self.metadata = metadata
    else:
        for x in metadata.columns:
            self.metadata[x] = pd.Series(metadata[x])

class Dandelion:
    def __init__(self, data=None, metadata=None, distance=None, edges=None, layout=None, graph=None):
        self.data = data
        self.metadata = metadata
        self.distance = distance
        self.edges = edges
        self.layout = layout
        self.graph = graph
        if self.data is not None:
            initialize_metadata(self)

    def copy(self):
        return copy.deepcopy(self)

    @staticmethod
    def isGZIP(filename):
        if filename.split('.')[-1] == 'gz':
            return True
        return False

    # Using HIGHEST_PROTOCOL is almost 2X faster and creates a file that
    # is ~10% smaller.  Load times go down by a factor of about 3X.
    def save(self, filename='dandelion_data.pkl'):
        if self.isGZIP(filename):
            f = gzip.open(filename, 'wb')
        else:
            f = open(filename, 'wb')
        pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.close()

    # Note that loading to a string with pickle.loads is about 10% faster
    # but probaly comsumes a lot more memory so we'll skip that for now.
    @classmethod
    def load(cls, filename='dandelion_data.pkl'):
        if cls.isGZIP(filename):
            f = gzip.open(filename, 'rb')
        else:
            f = open(filename, 'rb')
        n = pickle.load(f)
        f.close()
        return n

def convert_preprocessed_tcr_10x(file, prefix = None, save = None):
    """
    Parameters
    ----------
    file
        dandelion processed file
    prefix
        prefix to add to barcodes
    save
        save to specified location. Defaults to 'dandelion/data/'.
    """

    cr_annot = pd.read_csv(file, dtype = 'object')
    if prefix is not None:
        cr_annot['index']=[prefix+'_'+i for i in cr_annot['contig_id']]
    else:
        cr_annot['index']=[i for i in cr_annot['contig_id']]
    cr_annot.set_index('index', inplace = True)

    ddl_annot = pd.read_csv("{}/dandelion/data/tmp/{}_igblast.tsv".format(os.path.dirname(file), os.path.basename(file).split('_annotations.csv')[0]), sep = '\t', dtype = 'object')
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

    cellrangermap = {
        'sequence_id':'contig_id',
        'locus':'chain',
        'v_call_igblast':'v_gene',
        'd_call_igblast':'d_gene',
        'j_call_igblast':'j_gene',
        'productive':'productive',
        'junction_aa':'cdr3',
        'junction':'cdr3_nt'}

    for i in tqdm(cr_annot.index, desc = 'Matching and updating contig ids'):
        for key, value in cellrangermap.items():
            if cr_annot.loc[i, 'chain'] not in ['IGH', 'IGK', 'IGL', None]:
                cr_annot.loc[i, value] = ddl_annot.loc[i, key]
            else:
                cr_annot.loc[i, 'contig_id'] = ddl_annot.loc[i, 'sequence_id']

        if cr_annot.loc[i, 'cdr3'] is np.nan:
            cr_annot.loc[i, 'productive'] = None
        else:
            if cr_annot.loc[i, 'productive'] == 'T':
                cr_annot.loc[i, 'productive'] = 'True'
            else:
                cr_annot.loc[i, 'productive'] = 'False'

    cr_annot['barcode'] = [c.split('_contig')[0].split('-')[0] for c in cr_annot['contig_id']]
    if save is not None:
        cr_annot.to_csv(save, index = False, na_rep = 'None')
    else:
        cr_annot.to_csv("{}/dandelion/data/{}".format(os.path.dirname(file), os.path.basename(file)), index = False, na_rep = 'None')