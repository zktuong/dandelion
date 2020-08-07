#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 14:01:32
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-08-07 15:29:37

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
from changeo.IO import readGermlines
import warnings
try:
    from scanpy import logging as logg
except ImportError:
    pass

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
    """
    Reads in or copy dataframe object and set sequence_id as index without dropping.

    Parameters
    ----------
    obj : DataFrame, str
        file path to .tsv file or pandas DataFrame object.
    Returns
    -------
        pandas DataFrame object.
    """
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

def setup_metadata_(data):
    """
    A Dandelion class subfunction to initialize the `.metadata` slot.
    Parameters
    ----------
    data : DataFrame
        pandas DataFrame object.
    sep : tuple[int, str]
        A tuple containing how the clone groups should be extracted. None defaults to (0, '_').
    Returns
    -------
        pandas DataFrame object.
    """
    dat_h = data[data['locus'] == 'IGH']
    dict_h = dict(zip(dat_h['sequence_id'], dat_h['cell_id']))
    metadata_ = pd.DataFrame.from_dict(dict_h, orient = 'index', columns = ['cell_id'])
    metadata_.set_index('cell_id', inplace = True)    
    return(metadata_)

def setup_metadata(data, sep, clone_key = None):
    """
    A Dandelion class subfunction to initialize the `.metadata` slot.
    Parameters
    ----------
    data : DataFrame
        pandas DataFrame object.
    sep : tuple[int, str]
        A tuple containing how the clone groups should be extracted. None defaults to (0, '_').
    Returns
    -------
        pandas DataFrame object.
    """
    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    dat_h = data[data['locus'] == 'IGH']
    dat_l = data[data['locus'].isin(['IGK', 'IGL'])]
    clone_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h[clonekey])))
    clone_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l[clonekey])))
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
    metadata_[clonekey] = pd.Series(clones_list)
    metadata_ = metadata_[[clonekey]]
    if sep is None:
        scb = (0, '_')
    else:
        scb = (sep[0], sep[1])
    group = []
    # check if contain the separator
    x = list(metadata_[clonekey].str.contains(scb[1]))
    cl_ = list(metadata_[clonekey])
    for c in range(0, len(metadata_[clonekey])):
        if not x[c]:
            warnings.warn(UserWarning("\n\nSome/all clones do not contain '{}' as separator. \n".format(scb[1])))
            group.append(cl_[c])
        else:
            group.append(cl_[c].split(scb[1])[scb[0]])
    groupseries = dict(zip(metadata_.index, group))
    metadata_[str(clonekey)+'_group'] = pd.Series(groupseries)

    return(metadata_)

def retrieve_metadata(data, retrieve_id, split_heavy_light, collapse):
    """
    A Dandelion class subfunction to populate the `.metadata` slot.
    Parameters
    ----------
    data : DataFrame
        pandas DataFrame object.
    retrieve_id : str
        column name in `.data` slot.
    split_heavy_light : bool
        Returns the retrieval splitted into two column for heavy and light. Default is True. False combines the retrieval into a single column.
    collapse : bool
        Whether or not to collapse unique elements if duplicated. For example, different contigs and same sample id would then benefit from this option being set to True.
    Returns
    -------
        A dictionary with keys as cell_ids and records as retrieved value.
    """
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
    if not split_heavy_light:
        retrieval_list = {}
        for x in sub_metadata.index:
            if collapse:
                r_l = list(set(list(sub_metadata.loc[x, :])))
            else:
                r_l = list(sub_metadata.loc[x, :])
            r_l = sorted([y for y in r_l if str(y) != 'nan'])
            if len(r_l) > 1:
                r_l = r_l[1:]
            retrieval_list[x] = '|'.join(r_l)
        return(retrieval_list)
    else:
        heavy_retrieval_list = dict(sub_metadata['heavy'])
        for k, r in heavy_retrieval_list.items():
            if isinstance(r, pd.Series):
                heavy_retrieval_list[k] = ','.join([str(x) for x in flatten(r.to_list())])
        light_retrieval_list = {}
        sub_metadata2 = sub_metadata.drop('heavy', axis = 1).drop_duplicates() # TODO: need to work out if i need to adjust this when there is a situation that it can't remove non duplicates?
        for x in sub_metadata2.index:
            if collapse:
                r_l = list(set(list(sub_metadata2.loc[x, :])))
            else:
                r_l = list(sub_metadata2.loc[x, :])
            r_l = sorted([y for y in r_l if str(y) != 'nan'])
            light_retrieval_list[x] = ['|'.join(x) if len(x) > 0 else np.nan for x in [r_l]][0]
        return(heavy_retrieval_list, light_retrieval_list)

def initialize_metadata(self, retrieve = None, isotype_dict = None, split_heavy_light = True, collapse = False, clones_sep = None, clone_key = None):
    """
    A Dandelion function to update and populate the `.metadata` slot.

    Parameters
    ----------
    self : Dandelion
        Dandelion object
    retrieve : str
        column name in `.data` slot.
    split_heavy_light : bool
        Returns the retrieval splitted into two column for heavy and light. Default is True. False combines the retrieval into a single column.
    collapse : bool
        Whether or not to collapse unique elements if duplicated. For example, different contigs and same sample id would then benefit from this option being set to True.
    clones_sep : tuple[int, str]
        A tuple containing how the clone groups should be extracted. None defaults to (0, '_').
    Returns
    -------
        Dandelion object with `.metadata` slot initialized.
    """
    dat = load_data(self.data)
    for x in ['cell_id', 'locus', 'c_call', 'umi_count']:
        if x not in dat.columns:
            raise KeyError ("Please check your object. %s is not in the columns of input data." % x)
    
    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    metadata_status = self.metadata
    if metadata_status is None:
        if clonekey in dat.columns:
            self.metadata = setup_metadata(dat, clones_sep, clonekey)
        else:
            self.metadata = setup_metadata_(dat)
    else:
        self.metadata = self.metadata

    if 'sample_id' in dat.columns:
        samp_id = retrieve_metadata(dat, 'sample_id', False, True)
    if 'v_call_genotyped' in dat.columns:
        heavy_v_call, light_v_call = retrieve_metadata(dat, 'v_call_genotyped', True, False)
    else:
        heavy_v_call, light_v_call = retrieve_metadata(dat, 'v_call', True, False)
    heavy_j_call, light_j_call = retrieve_metadata(dat, 'j_call', True, False)
    heavy_c_call, light_c_call = retrieve_metadata(dat, 'c_call', True, False)
    heavy_umi, light_umi = retrieve_metadata(dat, 'umi_count', True, False)
    heavy_status, light_status = retrieve_metadata(dat, 'locus', True, False)
    status = pd.DataFrame([heavy_status, light_status], index = ['heavy', 'light']).T
    for i in status.index:
        try:
            status.at[i, 'status'] = status.loc[i,'heavy']+' + '+status.loc[i,'light']
        except:
            status.at[i, 'status'] = status.loc[i,'heavy'] + '_only'
    if isotype_dict is None:
        conversion_dict = {'igha1':'IgA', 'igha2':'IgA', 'ighm':'IgM', 'ighd':'IgD', 'ighe':'IgE', 'ighg1':'IgG', 'ighg2':'IgG', 'ighg3':'IgG', 'ighg4':'IgG', 'igkc':'IgK', 'iglc1':'IgL', 'iglc2':'IgL', 'iglc3':'IgL', 'iglc4':'IgL', 'iglc5':'IgL', 'iglc6':'IgL', 'iglc7':'IgL', 'igha':'IgA', 'ighg':'IgG', 'iglc':'IgL', 'nan':np.nan, np.nan:np.nan} # the key for IgG being igh is on purpose because of how the counter works
    else:
        conversion_dict = isotype_dict
    isotype = {}
    for k in heavy_c_call:
        if heavy_c_call[k] == heavy_c_call[k]:
            if ',' in heavy_c_call[k]:
                # iso_d = defaultdict(int)
                # for c in heavy_c_call[k].lower():
                #     iso_d[c] += 1
                isotype[k] = ','.join([str(z) for z in [conversion_dict[y.lower()] for y in set([re.sub('[0-9]', '', x) for x in heavy_c_call[k].split(',')])]])
            else:
                isotype[k] = conversion_dict[heavy_c_call[k].lower()]
        else:
            isotype[k] = heavy_c_call[k]
    lightchain = {}
    for k in light_c_call:
        if light_c_call[k] == light_c_call[k]:
            if '|' in light_c_call[k]:
                if ',' in light_c_call[k]:
                    lc_y = []
                    for x in light_c_call[k].lower().split('|'):
                        iso_d = defaultdict(int)
                        for c in x:
                            iso_d[c] += 1
                        lc_y.append(re.sub(',|[0-9]', '', ''.join([k_ for k_,v_ in iso_d.items() if v_ >= 2])))
                    lightchain[k] = [conversion_dict[y] for y in lc_y]
                else:
                    lightchain[k] = '|'.join([conversion_dict[x] for x in light_c_call[k].lower().split('|')])
            else:
                if ',' in light_c_call[k]:
                    iso_d = defaultdict(int)
                    for c in light_c_call[k].lower():
                        iso_d[c] += 1
                    lightchain[k] = conversion_dict[re.sub(',|[0-9]', '', ''.join([k_ for k_,v_ in iso_d.items() if v_ >= 2]))]
                else:
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
    productive = retrieve_metadata(dat, 'productive', False, False)
    if 'sample_id' in dat.columns:
        self.metadata['sample_id'] = pd.Series(samp_id)
    self.metadata['isotype'] = pd.Series(isotype)
    self.metadata['lightchain'] = pd.Series(lightchain)
    self.metadata['status'] = pd.Series(status['status'])
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
    self.metadata['vdj_status'] = pd.Series(multi)
    # return this in this order
    if metadata_status is None:
        if clonekey in self.data.columns:
            if 'sample_id' in self.data.columns:
                self.metadata = self.metadata[['sample_id', str(clonekey), str(clonekey)+'_group', 'isotype', 'lightchain', 'status', 'vdj_status', 'productive',  'umi_counts_heavy', 'umi_counts_light', 'v_call_heavy', 'v_call_light', 'j_call_heavy', 'j_call_light', 'c_call_heavy', 'c_call_light']]
            else:
                self.metadata = self.metadata[[str(clonekey), str(clonekey)+'_group', 'isotype', 'lightchain', 'productive', 'status', 'vdj_status', 'umi_counts_heavy', 'umi_counts_light',  'v_call_heavy', 'v_call_light', 'j_call_heavy', 'j_call_light', 'c_call_heavy', 'c_call_light']]
        else:
            if 'sample_id' in self.data.columns:
                self.metadata = self.metadata[['sample_id', 'isotype', 'lightchain', 'status', 'vdj_status', 'productive',  'umi_counts_heavy', 'umi_counts_light', 'v_call_heavy', 'v_call_light', 'j_call_heavy', 'j_call_light', 'c_call_heavy', 'c_call_light']]
            else:
                self.metadata = self.metadata[['isotype', 'lightchain', 'productive', 'status', 'vdj_status', 'umi_counts_heavy', 'umi_counts_light',  'v_call_heavy', 'v_call_light', 'j_call_heavy', 'j_call_light', 'c_call_heavy', 'c_call_light']]

    # new function to retrieve non-standard columns
    if retrieve is not None:
        if retrieve in dat.columns:
            if not split_heavy_light:
                if collapse:
                    retrieve_dict = retrieve_metadata(dat, retrieve, False, True)
                else:
                    retrieve_dict = retrieve_metadata(dat, retrieve, False, False)
                self.metadata[str(retrieve)] = pd.Series(retrieve_dict)
            else:
                if collapse:
                    h_retrieve_dict, l_retrieve_dict = retrieve_metadata(dat, retrieve, True, True)
                else:
                    h_retrieve_dict, l_retrieve_dict = retrieve_metadata(dat, retrieve, True, False)
                self.metadata[str(retrieve)+'_heavy'] = pd.Series(h_retrieve_dict)
                self.metadata[str(retrieve)+'_light'] = pd.Series(l_retrieve_dict)
        else:
            raise KeyError('Unknown column : \'%s\' to retrieve.' % retrieve)

class Dandelion:
    def __init__(self, data=None, metadata=None, germline = None, distance=None, edges=None, layout=None, graph=None, initialize = True, **kwargs):
        self.data = data
        self.metadata = metadata        
        self.distance = distance
        self.edges = edges
        self.layout = layout
        self.graph = graph
        self.threshold = None
        self.germline = {}
        if germline is not None:
            self.germline.update(germline)

        if os.path.isfile(str(self.data)):
            self.data = load_data(self.data)

        if self.data is not None:
            if initialize is True:
                initialize_metadata(self, **kwargs) 
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
    
    def update_germline(self, corrected = None, germline = None, org = 'human'):
        """
        Update germline reference with corrected sequences and store in Dandelion object.

        Parameters
        ----------
        self : Dandelion
            Dandelion object.
        corrected : dict, str, optional
            dictionary of corrected germline sequences or file path to corrected germline sequences fasta file.
        germline : str, optional
            path to germline database folder. Defaults to `$GERMLINE` environmental variable.
        org : str
            organism of reference folder. Default is 'human'.
        Returns
        -------
            updated germline reference diciontary in `.germline` slot.
        """
        start = logg.info('Updating germline reference')
        env = os.environ.copy()
        if germline is None:
            try:
                gml = env['GERMLINE']
            except:
                raise OSError('Environmental variable GERMLINE must be set. Otherwise, please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files.')
            gml = gml+'imgt/'+org+'/vdj/'
        else:
            if os.path.isdir(germline):
                gml = germline
            elif type(germline) is not list:
                germline_ = [germline]
                if len(germline_) < 3:
                    raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                else:
                    gml = []
                    for x in germline_:
                        if not x.endswith('.fasta'):
                            raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                        gml.append(x)
            elif type(germline) is list:
                if len(germline) < 3:
                    raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                else:
                    gml = []
                    for x in germline:
                        if not x.endswith('.fasta'):
                            raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                        gml.append(x)
            
        if type(gml) is not list:
            gml = [gml]

        germline_ref = readGermlines(gml)
        if corrected is not None:
            if type(corrected) is dict:
                personalized_ref_dict = corrected
            elif os.path.isfile(str(corrected)):
                personalized_ref_dict = readGermlines([corrected])
            # update with the personalized germline database
            germline_ref.update(personalized_ref_dict)
        else:
            pass

        self.germline.update(germline_ref)
        logg.info(' finished', time=start,
        deep=('Updated Dandelion object: \n'
        '   \'germline\', updated germline reference\n'))

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
    file : str
        file path to .tsv file.
    prefix : str
        prefix to add to barcodes. Ignored if left as None.
    save : str
        file path to save location. Defaults to 'dandelion/data/' if left as None.
    """

    cr_annot = load_data(file)
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
