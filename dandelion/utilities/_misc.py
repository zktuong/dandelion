#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 14:01:32
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-06-08 12:00:11

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
import warnings

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

def setup_metadata(data, sep):
        dat_h = data[data['locus'] == 'IGH']
        dat_l = data[data['locus'].isin(['IGK', 'IGL'])]        
        clone_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h['clone_id'])))
        clone_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l['clone_id'])))    
        metadata_ = pd.DataFrame.from_dict(clone_h, orient = 'index', columns = ['cell_id', 'heavy'])
        metadata_.set_index('cell_id', inplace = True)
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
        metadata_['light'] = pd.Series(light_clone_tree2)
        tmp_dat = metadata_['light'].apply(pd.Series)
        tmp_dat.columns = ['light_' + str(c) for c in tmp_dat.columns]
        metadata_ = metadata_.merge(tmp_dat, left_index = True, right_index = True)
        metadata_ = metadata_[['heavy'] + [str(c) for c in tmp_dat.columns]]
        clones_list = {}
        for x in metadata_.index:
            cl = list(set(list(metadata_.loc[x, :])))
            cl = sorted([y for y in cl if str(y) != 'nan'])
            if len(cl) > 1:
                cl = cl[1:]
            clones_list[x] = ','.join(cl)
        metadata_['clone_id'] = pd.Series(clones_list)
        metadata_ = metadata_[['clone_id']]
        if sep is None:
            scb = (0, '_')
        else:
            scb = (sep[0], sep[1])
        group = []
        
        # check if contain the separator
        x = list(metadata_['clone_id'].str.contains(scb[1]))
        cl_ = list(metadata_['clone_id'])
        for c in range(0, len(metadata_['clone_id'])):
            if not x[c]:
                warnings.warn(UserWarning("\n\nSome/all clones do not contain '{}' as separator. \n".format(scb[1])))
                group.append(cl_[c])
            else:
                group.append(cl_[c].split(scb[1])[scb[0]])
        groupseries = dict(zip(metadata_.index, group))
        metadata_['clone_group_id'] = pd.Series(groupseries)
        return(metadata_)

def retrieve_metadata(data, retrieve_id, combined):
        dat_h = data[data['locus'] == 'IGH']
        dat_l = data[data['locus'].isin(['IGK', 'IGL'])]
        
        retrieve_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h[retrieve_id])))
        retrieve_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l[retrieve_id])))
        
        sub_metadata = pd.DataFrame.from_dict(retrieve_h, orient = 'index', columns = ['cell_id', 'heavy'])
        sub_metadata.set_index('cell_id', inplace = True)

        light_retrieval_tree = Tree()
        for key, value in retrieve_l.items():
            k, v = value
            light_retrieval_tree[k][key] = v
        light_retrieval_tree2 = Tree()
        for g in light_retrieval_tree:
            second_key = []
            for k2 in light_retrieval_tree[g].keys():
                second_key.append(k2)
            second_key = list(set(second_key))
            second_key_dict = dict(zip(second_key, range(0,len(second_key))))
            for key, value in light_retrieval_tree[g].items():
                light_retrieval_tree2[g][second_key_dict[key]] = value
        sub_metadata['light'] = pd.Series(light_retrieval_tree2)
        tmp_dat = sub_metadata['light'].apply(pd.Series)
        tmp_dat.columns = ['light_' + str(c) for c in tmp_dat.columns]
        sub_metadata = sub_metadata.merge(tmp_dat, left_index = True, right_index = True)
        sub_metadata = sub_metadata[['heavy'] + [str(c) for c in tmp_dat.columns]]
        if combined:
            retrieval_list = {}
            for x in sub_metadata.index:
                r_l = list(set(list(sub_metadata.loc[x, :])))
                r_l = sorted([y for y in r_l if str(y) != 'nan'])
                if len(r_l) > 1:
                    r_l = r_l[1:]
                retrieval_list[x] = ','.join(r_l)
            return(retrieval_list)
        else:
            heavy_retrieval_list = dict(sub_metadata['heavy'])
            light_retrieval_list = {}
            sub_metadata2 = sub_metadata.drop('heavy', axis = 1)
            for x in sub_metadata2.index:
                r_l = list(set(list(sub_metadata2.loc[x, :])))
                r_l = sorted([y for y in r_l if str(y) != 'nan'])
                light_retrieval_list[x] = [','.join(x) if len(x) > 0 else np.nan for x in [r_l]][0]
            return(heavy_retrieval_list, light_retrieval_list)

def initialize_metadata(self, retrieve = None, collapse = False, clones_sep = None):
    # a quick way to retrieve the 'meta data' for each cell that transfer into obs lot in scanpy later
    
    dat = load_data(self.data)
    for x in ['cell_id', 'locus', 'c_call', 'umi_count']:
        if x not in dat.columns:
            raise KeyError ("Please check your object. %s is not in the columns of input data." % x)
    
    if 'clone_id' in dat. columns:
        self.metadata = setup_metadata(dat, clones_sep)

        if 'sample_id' in dat.columns:
            samp_id = retrieve_metadata(dat, 'sample_id', True)

        if 'v_call_genotyped' in dat.columns:
            heavy_v_call, light_v_call = retrieve_metadata(dat, 'v_call_genotyped', False)
        else:
            heavy_v_call, light_v_call = retrieve_metadata(dat, 'v_call', False)
        heavy_j_call, light_j_call = retrieve_metadata(dat, 'j_call', False)
        heavy_c_call, light_c_call = retrieve_metadata(dat, 'c_call', False)
        heavy_umi, light_umi = retrieve_metadata(dat, 'umi_count', False)
    
        conversion_dict = {'igha1':'IgA', 'igha2':'IgA', 'ighm':'IgM', 'ighd':'IgD', 'ighe':'IgE', 'ighg1':'IgG', 'ighg2':'IgG', 'ighg3':'IgG', 'ighg4':'IgG', 'igkc':'IgK', 'iglc1':'IgL', 'iglc2':'IgL', 'iglc3':'IgL', 'iglc4':'IgL', 'iglc5':'IgL', 'iglc6':'IgL', 'iglc7':'IgL'}
        isotype = {}
        for k in heavy_c_call:
            if heavy_c_call[k] == heavy_c_call[k]:
                isotype[k] = conversion_dict[heavy_c_call[k].lower()]
            else:
                isotype[k] = heavy_c_call[k]
        lightchain = {}
        for k in light_c_call:
            if light_c_call[k] == light_c_call[k]:
                lightchain[k] = conversion_dict[light_c_call[k].lower()]
            else:
                lightchain[k] = light_c_call[k]
        
        for k in heavy_v_call:
            heavy_v_call[k] = ''.join([','.join(list(set([re.sub('[*][0-9][0-9]', '', str(heavy_v_call[k]))][0].split(','))))])
        for k in heavy_j_call:
            heavy_j_call[k] = ''.join([','.join(list(set([re.sub('[*][0-9][0-9]', '', str(heavy_j_call[k]))][0].split(','))))])
        for k in light_v_call:
            light_v_call[k] = ''.join([','.join(list(set([re.sub('[*][0-9][0-9]', '', str(light_v_call[k]))][0].split(','))))])
        for k in light_j_call:
            light_j_call[k] = ''.join([','.join(list(set([re.sub('[*][0-9][0-9]', '', str(light_j_call[k]))][0].split(','))))])
    
        productive = retrieve_metadata(dat, 'productive', True)
        if 'sample_id' in dat.columns:
            self.metadata['sample_id'] = pd.Series(samp_id)
        self.metadata['isotype'] = pd.Series(isotype)
        self.metadata['lightchain'] = pd.Series(lightchain)
        
        self.metadata['productive'] = pd.Series(productive)
        self.metadata['umi_counts_heavy'] = pd.Series(heavy_umi)
        self.metadata['umi_counts_light'] = pd.Series(light_umi)
        
        self.metadata['c_call_heavy'] = pd.Series(heavy_c_call)
        self.metadata['c_call_light'] = pd.Series(light_c_call)
        self.metadata['v_call_heavy'] = pd.Series(heavy_v_call)
        self.metadata['v_call_light'] = pd.Series(light_v_call)
        self.metadata['j_call_heavy'] = pd.Series(heavy_j_call)
        self.metadata['j_call_light'] = pd.Series(light_j_call)
    
        multi = {}
        for i in self.metadata.index:        
            try:
                hv_ = self.metadata.loc[i, 'v_call_heavy'].split(',')
            except:
                hv_ = self.metadata.loc[i, 'v_call_heavy']
            try:
                hj_ = self.metadata.loc[i, 'j_call_heavy'].split(',')
            except:
                hj_ = self.metadata.loc[i, 'j_call_heavy']
            try:
                lv_ = self.metadata.loc[i, 'v_call_light'].split(',')
            except:
                lv_ = self.metadata.loc[i, 'v_call_light']
            try:            
                lj_ = self.metadata.loc[i, 'v_call_light'].split(',')
            except:
                lv_ = self.metadata.loc[i, 'v_call_light']
    
            multi_ = []
            if len(hv_) > 1:
                multi_.append(['Multi_heavy_v'])
            if len(hj_) > 1:
                multi_.append(['Multi_heavy_j'])
            if len(lv_) > 1:
                multi_.append(['Multi_light_v'])
            if len(lj_) > 1:
                multi_.append(['Multi_light_j'])
            if len(multi_) < 1:
                multi_.append(['Single'])
    
            multi[i] = ','.join(list(flatten(multi_)))
        self.metadata['multi_status'] = pd.Series(multi)
        
        # return this in this order
        if 'sample_id' in dat.columns:
            self.metadata = self.metadata[['sample_id', 'clone_id', 'clone_group_id', 'isotype', 'lightchain', 'productive', 'umi_counts_heavy', 'umi_counts_light', 'multi_status', 'v_call_heavy', 'v_call_light', 'j_call_heavy', 'j_call_light', 'c_call_heavy', 'c_call_light']]
        else:
            self.metadata = self.metadata[['clone_id', 'clone_group_id', 'isotype', 'lightchain', 'productive', 'umi_counts_heavy', 'umi_counts_light', 'multi_status', 'v_call_heavy', 'v_call_light', 'j_call_heavy', 'j_call_light', 'c_call_heavy', 'c_call_light']]
        # new function to retrieve non-standard columns
        if retrieve is not None:
            if retrieve in dat.columns:
                if collapse:
                    retrieve_dict = retrieve_metadata(dat, retrieve, True)
                    self.metadata[str(retrieve)] = pd.Series(retrieve_dict)
                else:
                    h_retrieve_dict, l_retrieve_dict = retrieve_metadata(dat, retrieve, False)
                    self.metadata[str(retrieve)+'_heavy'] = pd.Series(h_retrieve_dict)
                    self.metadata[str(retrieve)+'_light'] = pd.Series(l_retrieve_dict)
            else:
                raise KeyError('Unknown column : \'%s\' to retrieve.' % retrieve)
    else:
        pass

class Dandelion:
    def __init__(self, data=None, metadata=None, distance=None, edges=None, layout=None, graph=None):
        self.data = data
        self.metadata = metadata
        self.distance = distance
        self.edges = edges
        self.layout = layout
        self.graph = graph
        self.threshold = None
        
        if os.path.isfile(str(self.data)):
            self.data = load_data(self.data)
        
        if self.data is not None:
            initialize_metadata(self)
            self.n_contigs = self.data.shape[0]
            try:
                self.n_obs = self.metadata.shape[0]
            except:
                self.n_obs = 0
    
    def _gen_repr(self, n_obs, n_contigs) -> str:
        descr = f"Dandelion class object with n_obs = {n_obs} and n_contigs = {n_contigs}"
        for attr in ["data", "metadata", "distance", "edges"]:
            try:
                keys = getattr(self, attr).keys()
            except:
                keys = []
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(list(keys))[1:-1]}"
            else:
                descr += f"\n    {attr}: {str(None)}"
        descr += f"\n    layout: {str(self.layout).strip('<>')}"
        descr += f"\n    graph: {str(type(self.graph)).strip('<>')}"
        descr += f"\n    threshold: {self.threshold}"
        return descr

    def __repr__(self) -> str:
        return self._gen_repr(self.n_obs, self.n_contigs)

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