#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-05-13 23:22:18
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-08-07 23:25:15

import os
import sys
import scanpy as sc
import pandas as pd
from pandas import DataFrame
import numpy as np
from tqdm import tqdm
from ..utilities._utilities import *
from collections import defaultdict
from itertools import groupby
from scipy.spatial.distance import pdist, squareform
from distance import hamming
import re
import math
import scipy
from scipy.sparse.csgraph import minimum_spanning_tree
import Levenshtein
import networkx as nx
import igraph
from time import sleep
import copy
import functools
try:
    from scanpy import logging as logg
except ImportError:
    pass
import warnings
from subprocess import run
import multiprocessing
from changeo.Gene import getGene

def find_clones(self, identity=0.85, clustering_by = None, by_alleles = None, write_out = False, outdir=None, outFilePrefix=None, key_added = None):
    """
    Find clones based on heavy chain and light chain CDR3 junction hamming distance.

    Parameters
    ----------
    self : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones have been determined.
    identity : float
        Junction similarity parameter. Default 0.85
    clustering_by : str, optional
        modes for clustering: 'nt' or 'aa'. None defaults to 'aa'.
    by_alleles : bool, optional
        Whether or not to collapse alleles to genes. None defaults to True.
    write_out : bool
        If True, will write out airr/changeo file with clone_id column (default is False). file path and file name is determined by outdir and outFilePrefix options..
    outdir : str, optional
        If specified, outfile will be in this location. None defaults to 'dandelion/data'.
    outFilePrefix : str, optional
        If specified, the outfile name will have this prefix. None defaults to 'dandelion_find_clones'
    key_added : str, optional
        If specified, this will be the column name for clones. None defaults to 'clone_id'
    Returns
    -------
        `Dandelion` object with clone_id annotated in `.data` slot and `.metadata` initialized.
    """
    start = logg.info('Finding clones')
    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    else:
        dat = load_data(self)
    dat_heavy = dat[dat['locus'] == 'IGH']
    pd.set_option('mode.chained_assignment', None)

    if key_added is None:
        clone_key = 'clone_id'
    else:
        clone_key = key_added

    # retrieve the J genes and J genes
    if by_alleles is None or ~by_alleles:
        if 'v_call_genotyped' in dat_heavy.columns:
            V = [re.sub('[*][0-9][0-9]', '', v) for v in dat_heavy['v_call_genotyped']]
        else:
            V = [re.sub('[*][0-9][0-9]', '', v) for v in dat_heavy['v_call']]
        J = [re.sub('[*][0-9][0-9]', '', j) for j in dat_heavy['j_call']]
    elif by_alleles:
        if 'v_call_genotyped' in dat_heavy.columns:
            V = [v for v in dat_heavy['v_call_genotyped']]
        else:
            V = [v for v in dat_heavy['v_call']]
        J = [j for j in dat_heavy['j_call']]
    else:
        raise ValueError("by_alleles only accepts boolean values or None, with None defaulting to False.")

    # collapse the alleles to just genes
    V = [','.join(list(set(v.split(',')))) for v in V]
    J = [','.join(list(set(j.split(',')))) for j in J]

    if clustering_by is None or clustering_by is 'aa':
        junction = dict(zip(dat_heavy.index, dat_heavy['junction_aa']))
        if 'junction_aa_length' not in dat_heavy.columns:
            junction_length = [len(str(l)) for l in dat_heavy['junction_aa']]
        else:
            junction_length = [l for l in dat_heavy['junction_aa_length']]
        junction_length_dict = dict(zip(dat_heavy.index, junction_length))
    elif clustering_by == 'nt':
        junction = dict(zip(dat_heavy.index, dat_heavy['junction']))
        if 'junction_length' not in dat_heavy.columns:
            junction_length = [len(str(l)) for l in dat_heavy['junction']]
        else:
            junction_length = [l for l in dat_heavy['junction_length']]
        junction_length_dict = dict(zip(dat_heavy.index, junction_length))
    else:
        raise ValueError("clustering_by only accepts string values 'aa', 'nt' or None, with None defaulting to 'aa'.")
    # Create a dictionary and group sequence ids with same V and J genes
    V_J = dict(zip(dat_heavy.index, zip(V,J)))
    vj_grp = defaultdict(list)
    for key, val in sorted(V_J.items()):
        vj_grp[val].append(key)
    # and now we split the groups based on lengths of the junctions
    vj_len_grp = Tree()
    junction_grp = Tree()
    for g in vj_grp:
        # first, obtain what's the unique lengths
        jlen = []
        for contig_id in vj_grp[g]:
            jlen.append(junction_length_dict[contig_id])
        setjlen = list(set(jlen))
        # then for each unique length, we add the contigs to a new tree if it matches the length
        # and also the actual junction sequences into a another one
        for s in setjlen:
            for contig_id in vj_grp[g]:
                jlen_ = junction_length_dict[contig_id]
                if jlen_ == s:
                    vj_len_grp[g][s][contig_id].value = 1
                    junction_grp[g][s][junction[contig_id]].value = 1
                    for c in [contig_id]:
                        vj_len_grp[g][s][c] = junction[c]
    clones = Tree()
    # for each junction group, calculate the hamming distance matrix
    for g in tqdm(junction_grp, desc = 'Finding clones based on heavy chains '):
        for l in junction_grp[g]:
            junction_ = list(junction_grp[g][l])
            tdarray = np.array(junction_).reshape(-1,1)
            d_mat = squareform(pdist(tdarray, lambda x,y: hamming(x[0],y[0])))
            # then calculate what the acceptable threshold is for each length of sequence
            tr = math.floor(int(l)*(1-identity))
            # convert diagonal and upper triangle to zeroes
            d_mat = np.tril(d_mat)
            np.fill_diagonal(d_mat, 0)
            # get the coordinates/indices of junctions to match against the threshold later
            indices_temp = []
            indices = []
            indices_temp = [list(x) for x in np.tril_indices_from(d_mat)]
            indices = list(zip(indices_temp[0], indices_temp[1]))
            # if there's more than 1 contig, remove the diagonal
            if len(indices) > 1:
                for pairs in indices:
                    # remove diagonals
                    if pairs[0] == pairs[1]:
                        indices.remove(pairs)
            indices_j = []
            # use the coordinates/indices to retrieve the junction sequences
            for p in range(0, len(indices)):
                a1, b1 = indices[p]
                indices_j.append(junction_[a1])
                indices_j.append(junction_[b1])
            # retain only the unique sequences
            indices_j_f = list(set(indices_j))
            # convert the distance matrix to coordinate (source) and distance (target) and create it as a dictionary
            source, target = d_mat.nonzero()
            source_target = list(zip(source.tolist(), target.tolist()))
            if len(source) == 0 & len(target) == 0:
                source_target = list([(0,0)])
            dist = {}
            for st in source_target:
                dist.update({st:d_mat[st]})
            cm1 = []
            cm2 = []
            cm3 = []
            # now to calculate which contigs to group
            tr2 = min(dist.values())
            if tr2 <= tr:
                for d in dist:
                    # if the distance is equal to the minimum distance, which is lesser than or equal to the threshold. I will add the sequence to cm1, meaning they passed and should be a clone.
                    if dist[d] == tr2:
                        cm1.append(d)
                    else:
                        # otherwise add into cm3, which later on will get split into individual clones, despite being part of the acceptable threshold for this group of sequences
                        cm3.append(d)
            else:
                # if the value is greater than the the minimum threshold, I will add the entire group into cm3, so they shouldn't be grouped together as a clone
                for d in dist:
                    cm3.append(d)
            # This is based on a simplest scenario e.g. in 3x3 distance matrix where the acceptable distance threshold is 2.
            # SeqA is 1 sequence different from SeqB
            # SeqB is 2 sequences different from SeqC
            # SeqC can only be 2 sequences different from SeqA otherwise it will violate the pair-wise distance matrix calculation
            # Because the distance between SeqA and SeqB is below the threshold, and is essentially the smallest value, they will be grouped in cm1
            # While SeqC is below the acceptable threshold, really it's more different than SeqA/B 'clone_id'. Hence i'm splitting up here.
            # If there's more than 3 sequences in a particular distance matrix, it gets slightly complicated.
            # so i repeat the calculation to only catch those with more than 3 sequences
            if len(dist) > 3:
                for d in dist:
                    # This time, i will rely on the main threshold value; if the distance lesser than the overall threshold, i will add it to cm2.
                    if dist[d] < tr:
                        cm2.append(d)
                    else:
                        cm3.append(d)
                # any sequences that were already binned above in cm1, it will be removed from here.
                cm2 = list(set(cm2) ^ set(cm1))
                cm3 = list(set(cm3) ^ set(cm1))

            # now to actually retrieve the junction/contigs properly through a series of appending and removing
            j_list1 = []
            if len(cm1) > 0:
                for i in range(0, len(cm1)):
                    a, b = cm1[i]
                    j_list1.append(junction_[a])
                    j_list1.append(junction_[b])
                j_list1 = list(set(j_list1))

            j_list2 = []
            if len(cm3) > 0:
                for i in range(0, len(cm2)):
                    a, b = cm2[i]
                    j_list2.append(junction_[a])
                    j_list2.append(junction_[b])
                j_list2 = list(set(j_list2))

            j_list3 = []
            if len(cm3) > 0:
                for i in range(0, len(cm3)):
                    a, b = cm3[i]
                    j_list3.append(junction_[a])
                    j_list3.append(junction_[b])
                j_list3 = list(set(j_list3))

            # this is to catch the overlapping junctions appearing in the lists because of the sequential appending
            for jl3_1 in j_list1:
                if jl3_1 in j_list3:
                    j_list3.remove(jl3_1)
            for jl2_1 in j_list1:
                if jl2_1 in j_list2:
                    j_list2.remove(jl2_1)
            for jl3_2 in j_list2:
                if jl3_2 in j_list3:
                    j_list3.remove(jl3_2)

            j_list3 = [i.split() for i in list(set(j_list3))]

            if len(j_list1) > 0:
                clones[g][l][str(0)] = j_list1
            if len(j_list2) > 0:
                clones[g][l][str(1)] = j_list2
            if len(j_list3) > 0:
                for c in range(0, len(j_list3)):
                    # the +2 here is so that the numbers come up after 1. It doesn't matter because i will reformat the clone ID numbers later
                    clones[g][l][str(c+2)] = j_list3[c]

    clone_dict = {}
    # now to retrieve the contig ids that are grouped together
    cid = Tree()
    for g in clones:
        for l in clones[g]:
            # retrieve the clone 'numbers'
            for c in clones[g][l]:
                grp_junction = clones[g][l][c]
                for key, value in vj_len_grp[g][l].items():
                    if value in grp_junction:
                        cid[g][l][c][key].value = 1
    # rename clone ids - get dictionaries step by step
    first_key = []
    for k1 in cid.keys():
        first_key.append(k1)
    first_key = list(set(first_key))
    first_key_dict = dict(zip(first_key, range(1,len(first_key)+1)))
    # and now for the middle key
    for g in cid:
        second_key = []
        for k2 in cid[g].keys():
            second_key.append(k2)
        second_key = list(set(second_key))
        second_key_dict = dict(zip(second_key, range(1,len(second_key)+1)))
        for l in cid[g]:
            # and now for the last key
            third_key = []
            for k3 in cid[g][l].keys():
                third_key.append(k3)
            third_key = list(set(third_key))
            third_key_dict = dict(zip(third_key, range(1,len(third_key)+1)))
            for key, value in dict(cid[g][l]).items():
                vL = []
                for v in value:
                    if type(v) is int:
                        break
                    # instead of converting to another tree, i will just make it a dictionary
                    clone_dict[v] = str(first_key_dict[g])+'_'+str(second_key_dict[l])+'_'+str(third_key_dict[key])
    # add it to the original dataframes
    dat_heavy[clone_key] = pd.Series(clone_dict)
    hclone = dict(zip(dat_heavy['cell_id'], dat_heavy[clone_key]))
    hlclone = dict(zip(dat['sequence_id'], [hclone[c] for c in dat['cell_id']]))
    dat[clone_key] = pd.Series(hlclone)
    # repeat this process for the light chains within each clone, but only for those with more than 1 light chains in a clone
    dat_light = dat[~(dat['locus'] == 'IGH')]
    # retrieve the J genes and J genes
    for c in tqdm(list(set(dat_light[clone_key])), desc = 'Refining clone assignment based on light chain pairing '):
        dat_light_c = dat_light[dat_light[clone_key] == c]
        if dat_light_c.shape[0] > 1:
            if by_alleles is None or ~by_alleles:
                if 'v_call_genotyped' in dat_light_c.columns:
                    Vlight = [re.sub('[*][0-9][0-9]', '', v) for v in dat_light_c['v_call_genotyped']]
                else:
                    Vlight = [re.sub('[*][0-9][0-9]', '', v) for v in dat_light_c['v_call']]
                Jlight = [re.sub('[*][0-9][0-9]', '', j) for j in dat_light_c['j_call']]
            elif by_alleles:
                if 'v_call_genotyped' in dat_light_c.columns:
                    Vlight = [v for v in dat_light_c['v_call_genotyped']]
                else:
                    Vlight = [v for v in dat_light_c['v_call']]
                Jlight = [j for j in dat_light_c['j_call']]
            else:
                raise ValueError("by_alleles only accepts boolean values or None, with None defaulting to False.")
            # collapse the alleles to just genes
            Vlight = [','.join(list(set(v.split(',')))) for v in Vlight]
            Jlight = [','.join(list(set(j.split(',')))) for j in Jlight]
            if clustering_by is None or clustering_by is 'aa':
                junction = dict(zip(dat_light_c.index, dat_light_c['junction_aa']))
                if 'junction_aa_length' not in dat_light.columns:
                    junction_length = [len(str(l)) for l in dat_light_c['junction_aa']]
                else:
                    junction_length = [l for l in dat_light_c['junction_aa_length']]
                junction_length_dict = dict(zip(dat_light_c.index, junction_length))
            elif clustering_by == 'nt':
                junction = dict(zip(dat_light_c.index, dat_light_c['junction']))
                if 'junction_length' not in dat_light_c.columns:
                    junction_length = [len(str(l)) for l in dat_light_c['junction']]
                else:
                    junction_length = [l for l in dat_light_c['junction_length']]
                junction_length_dict = dict(zip(dat_light_c.index, junction_length))
            else:
                raise ValueError("clustering_by only accepts string values 'aa', 'nt' or None, with None defaulting to 'aa'.")

            # Create a dictionary and group sequence ids with same V and J genes
            V_Jlight = dict(zip(dat_light_c.index, zip(Vlight,Jlight)))
            vj_lightgrp = defaultdict(list)
            for key, val in sorted(V_Jlight.items()):
                vj_lightgrp[val].append(key)
            # and now we split the groups based on lengths of the junctions
            vj_len_lightgrp = Tree()
            junction_lightgrp = Tree()
            for g in vj_lightgrp:
                # first, obtain what's the unique lengths
                jlen = []
                for contig_id in vj_lightgrp[g]:
                    jlen.append(junction_length_dict[contig_id])
                setjlen = list(set(jlen))
                # then for each unique length, we add the contigs to a new tree if it matches the length
                # and also the actual junction sequences into a another one
                for s in setjlen:
                    for contig_id in vj_lightgrp[g]:
                        jlen_ = junction_length_dict[contig_id]
                        if jlen_ == s:
                            vj_len_lightgrp[g][s][contig_id].value = 1
                            junction_lightgrp[g][s][junction[contig_id]].value = 1
                            for c in [contig_id]:
                                vj_len_lightgrp[g][s][c] = junction[c]

            clones_light = Tree()
            for g in junction_lightgrp:
                for l in junction_lightgrp[g]:
                    junction_ = list(junction_lightgrp[g][l])
                    tdarray = np.array(junction_).reshape(-1,1)
                    d_mat = squareform(pdist(tdarray, lambda x,y: hamming(x[0],y[0])))
                    tr = math.floor(int(l)*(1-identity))
                    d_mat = np.tril(d_mat)
                    np.fill_diagonal(d_mat, 0)
                    indices_temp = []
                    indices = []
                    indices_temp = [list(x) for x in np.tril_indices_from(d_mat)]
                    indices = list(zip(indices_temp[0], indices_temp[1]))
                    if len(indices) > 1:
                        for pairs in indices:
                            if pairs[0] == pairs[1]:
                                indices.remove(pairs)
                    indices_j = []
                    for p in range(0, len(indices)):
                        a1, b1 = indices[p]
                        indices_j.append(junction_[a1])
                        indices_j.append(junction_[b1])
                    indices_j_f = list(set(indices_j))
                    source, target = d_mat.nonzero()
                    source_target = list(zip(source.tolist(), target.tolist()))
                    if len(source) == 0 & len(target) == 0:
                        source_target = list([(0,0)])
                    dist = {}
                    for st in source_target:
                        dist.update({st:d_mat[st]})
                    cm1 = []
                    cm2 = []
                    cm3 = []
                    tr2 = min(dist.values())
                    if tr2 <= tr:
                        for d in dist:
                            if dist[d] == tr2:
                                cm1.append(d)
                            else:
                                cm3.append(d)
                    else:
                        for d in dist:
                            cm3.append(d)

                    if len(dist) > 3:
                        for d in dist:
                            if dist[d] < tr:
                                cm2.append(d)
                            else:
                                cm3.append(d)
                        cm2 = list(set(cm2) ^ set(cm1))
                        cm3 = list(set(cm3) ^ set(cm1))
                    j_list1 = []
                    if len(cm1) > 0:
                        for i in range(0, len(cm1)):
                            a, b = cm1[i]
                            j_list1.append(junction_[a])
                            j_list1.append(junction_[b])
                        j_list1 = list(set(j_list1))
                    j_list2 = []
                    if len(cm3) > 0:
                        for i in range(0, len(cm2)):
                            a, b = cm2[i]
                            j_list2.append(junction_[a])
                            j_list2.append(junction_[b])
                        j_list2 = list(set(j_list2))
                    j_list3 = []
                    if len(cm3) > 0:
                        for i in range(0, len(cm3)):
                            a, b = cm3[i]
                            j_list3.append(junction_[a])
                            j_list3.append(junction_[b])
                        j_list3 = list(set(j_list3))
                    for jl3_1 in j_list1:
                        if jl3_1 in j_list3:
                            j_list3.remove(jl3_1)
                    for jl2_1 in j_list1:
                        if jl2_1 in j_list2:
                            j_list2.remove(jl2_1)
                    for jl3_2 in j_list2:
                        if jl3_2 in j_list3:
                            j_list3.remove(jl3_2)
                    j_list3 = [i.split() for i in list(set(j_list3))]
                    if len(j_list1) > 0:
                        clones_light[g][l][str(0)] = j_list1
                    if len(j_list2) > 0:
                        clones_light[g][l][str(1)] = j_list2
                    if len(j_list3) > 0:
                        for c in range(0, len(j_list3)):
                            clones_light[g][l][str(c+2)] = j_list3[c]
            clone_dict_light = {}
            cid_light = Tree()
            for g in clones_light:
                for l in clones_light[g]:
                    for c in clones_light[g][l]:
                        grp_junction = clones_light[g][l][c]
                        for key, value in vj_len_lightgrp[g][l].items():
                            if value in grp_junction:
                                cid_light[g][l][c][key].value = 1
            first_key = []
            for k1 in cid_light.keys():
                first_key.append(k1)
            first_key = list(set(first_key))
            first_key_dict = dict(zip(first_key, range(1,len(first_key)+1)))
            for g in cid_light:
                second_key = []
                for k2 in cid_light[g].keys():
                    second_key.append(k2)
                second_key = list(set(second_key))
                second_key_dict = dict(zip(second_key, range(1,len(second_key)+1)))
                for l in cid_light[g]:
                    third_key = []
                    for k3 in cid_light[g][l].keys():
                        third_key.append(k3)
                    third_key = list(set(third_key))
                    third_key_dict = dict(zip(third_key, range(1,len(third_key)+1)))
                    for key, value in dict(cid_light[g][l]).items():
                        vL = []
                        for v in value:
                            if type(v) is int:
                                break
                            clone_dict_light[v] = str(first_key_dict[g])+'_'+str(second_key_dict[l])+'_'+str(third_key_dict[key])
            lclones = list(clone_dict_light.values())
            # will just update the main dat directly
            if len(list(set(lclones))) > 1:
                lclones_dict = dict(zip(sorted(list(set(lclones))), [str(x) for x in range(1,len(list(set(lclones)))+1)]))
                renamed_clone_dict_light = {}
                for key, value in clone_dict_light.items():
                    renamed_clone_dict_light[key] = lclones_dict[value]
                dat.at[renamed_clone_dict_light.keys(), clone_key] = dat.loc[renamed_clone_dict_light.keys(), clone_key] + '_' + pd.Series(renamed_clone_dict_light)

    if os.path.isfile(str(self)):
        dat.to_csv("{}/{}_clone.tsv".format(os.path.dirname(self), os.path.basename(self).split('.tsv')[0]), sep = '\t', index = False)
    else:
        if write_out:
            if outdir is None:
                outDir = 'dandelion/data'
            else:
                if outdir.endswith('/'):
                    outDir = str(outdir).strip('/')
                else:
                    outDir = str(outdir)

            if not os.path.exist(outDir):
                os.makedirs(outDir)
            if outFilePrefix is not None:
                dat.to_csv("{}/{}_clone.tsv".format(outDir, str(outFilePrefix)), sep = '\t', index = None)
            elif outFilePrefix is None:
                dat.to_csv("{}/{}_clone.tsv".format(outDir, 'dandelion_find_clones'), sep = '\t', index = None)
        else:
            pass
    sleep(0.5)
    logg.info(' finished', time=start,
        deep=('Updated Dandelion object: \n'
        '   \'data\', contig-indexed clone table\n'
        '   \'metadata\', cell-indexed clone table\n'))
    if self.__class__ == Dandelion:
        if self.germline is not None:
            germline_ = self.germline
        else:
            germline_ = None
        if self.distance is not None:
            dist_ = self.distance
        else:
            dist_ = None
        if self.edges is not None:
            edge_ = self.edges
        else:
            edge_ = None
        if self.layout is not None:
            layout_ = self.layout
        else:
            layout_ = None
        if self.graph is not None:
            graph_ = self.graph
        else:
            graph_ = None
        if self.threshold is not None:
            threshold_ = self.threshold
        else:
            threshold_ = None
        if ('clone_id' in self.data.columns) and (clone_key is not None):
            self.__init__(data = dat, germline = germline_, distance = dist_, edges = edge_, layout = layout_, graph = graph_, initialize = True, retrieve = clone_key, split_heavy_light = False)
        elif ('clone_id' not in self.data.columns) and (clone_key is not None):
            self.__init__(data = dat, germline = germline_, distance = dist_, edges = edge_, layout = layout_, graph = graph_, initialize = True, clone_key = clone_key, retrieve = clone_key, split_heavy_light = False)
        else:
            self.__init__(data = dat, germline = germline_, distance = dist_, edges = edge_, layout = layout_, graph = graph_, initialize = True, clone_key = clone_key)
        self.threshold = threshold_
    else:
        out = Dandelion(data = dat, clone_key = clone_key, retrieve = clone_key, split_heavy_light = False)
        return(out)

def generate_network(self, distance_mode='simple', aa_or_nt=None, clone_key = None, clones_sep = None, weights = None, layout_option = None, *args, **kwds):
    """
    Generates a levenshtein distance network based on gapped full length sequences for heavy and light chain(s).
    The distance matrices are then combined into a singular matrix where a minimum spanning tree will be constructed per clone group specified by separator in `clones_sep` option.

    Parameters
    ----------
    data : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones have been determined.
    distance_mode : str
        The mode of calculating joint distance matrix for heavy and light chains. Default is 'simple'. If 'simple', a simple sum operation will be used. If 'weighted', depending on whether `weights` option is provided, it will scale each layer to range of 0 to 1 to bring the multiple layers of data into a single analysis.
    aa_or_nt : str, optional
        Option accepts 'aa', 'nt' or None, with None defaulting to 'aa'. Determines whether amino acid or nucleotide sequences will be used for calculating distances.
    clones_sep: tuple[int, str]
        A tuple containing how the clone groups should be extracted. None defaults to (0, '_').
    weights : tuple, optional
        A tuple containing weights to scale each layer. default is None where each layer is scaled evenly i.e. 1/number of layers.
    layout_option : str, optional
        choice of layout algorithm. None defaults to fruchterman reingold layout.
    *args and **kwds
        passed to `igraph.graph.layout <https://igraph.org/python/doc/igraph.Graph-class.html>`__.
    Returns
    ----------
        `Dandelion` object with `.distance`, `.edges`, `.layout`, `.graph` initialized.
    """
    start = logg.info('Generating network')
    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    else:
        dat = load_data(self)
    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key
    if clonekey not in dat.columns:
        raise TypeError('Data does not contain clone information. Please run find_clones.')

    # re-initiate a Dandelion class object
    out = Dandelion(dat)

    # calculate distance
    dat_h = dat[dat['locus'] == 'IGH']
    dat_l = dat[dat['locus'].isin(['IGK', 'IGL'])]
    if aa_or_nt is None or aa_or_nt is 'aa':
        seq_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h['sequence_alignment_aa'])))
        seq_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l['sequence_alignment_aa'])))
    elif aa_or_nt == 'nt':
        seq_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h['sequence_alignment'])))
        seq_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l['sequence_alignment'])))
    else:
        raise ValueError("aa_or_nt only accepts string values 'aa', 'nt' or None, with None defaulting to 'aa'.")

    # So first, create a data frame to hold all possible (full) sequences split by heavy (only 1 possible) and light (multiple possible)
    dat_seq = pd.DataFrame.from_dict(seq_h, orient = 'index', columns = ['cell_id', 'heavy'])
    dat_seq.set_index('cell_id', inplace = True)
    light_seq_tree = Tree()
    for key, value in seq_l.items():
        k, v = value
        light_seq_tree[k][key] = v
    light_seq_tree2 = Tree()
    for g in light_seq_tree:
        second_key = []
        for k2 in light_seq_tree[g].keys():
            second_key.append(k2)
        second_key = list(set(second_key))
        second_key_dict = dict(zip(second_key, range(0,len(second_key))))
        for key, value in light_seq_tree[g].items():
            light_seq_tree2[g][second_key_dict[key]] = value
    dat_seq['light'] = pd.Series(light_seq_tree2)
    tmp_dat = dat_seq['light'].apply(pd.Series)
    tmp_dat.columns = ['light_' + str(c) for c in tmp_dat.columns]
    dat_seq = dat_seq.merge(tmp_dat, left_index = True, right_index = True)
    dat_seq = dat_seq[['heavy'] + [str(c) for c in tmp_dat.columns]]

    # calculate a distance matrix for all vs all and this can be referenced later on to extract the distance between the right pairs
    dmat = Tree()
    for x in tqdm(dat_seq.columns, desc = 'Calculating distances... '):
        seq_list = []
        seq_list = [y for y in dat_seq[x]]
        tdarray = np.array(seq_list).reshape(-1,1)
        d_mat = squareform(pdist(tdarray,lambda x,y: Levenshtein.distance(x[0],y[0])))
        dmat[x] = d_mat
    dist_mat_list = [dmat[x] for x in dmat if type(dmat[x]) is np.ndarray]

    n_ = len(dist_mat_list)
    if distance_mode == 'simple':
        total_dist = np.sum(dist_mat_list,axis=0)
    if distance_mode == 'weighted':
        weighted_matrix = []
        if weights is None:
            for w in range(0, n_):
                weighted_matrix.append(1/n_ * dist_mat_list[w])
            total_dist = sum(weighted_matrix)
        else:
            if len(weights) == n_:
                for w in range(0, n_):
                    weighted_matrix.append(weights[w] * dist_mat_list[w])
                total_dist = sum(weighted_matrix)
            else:
                raise IndexError('Length of provided weights should be %s.' % int(n_))

    # generate edge list
    tmp_totaldist = pd.DataFrame(total_dist, index = out.metadata.index, columns = out.metadata.index)
    tmp_clusterdist = Tree()
    for i in out.metadata.index:
        cx = out.metadata.loc[i, str(clonekey)+'_group']
        tmp_clusterdist[cx][i].value = 1
    tmp_clusterdist2 = {}
    for x in tmp_clusterdist:
        tmp_clusterdist2[x] = list(tmp_clusterdist[x])
    cluster_dist = {}
    for x in tmp_clusterdist2:
        dist_mat_ = tmp_totaldist.loc[tmp_clusterdist2[x], tmp_clusterdist2[x]]
        s1, s2 = dist_mat_.shape
        if s1 > 1 and s2 >1:
            cluster_dist[x] = dist_mat_
    # to improve the visulisation and plotting efficiency, i will build a minimum spanning tree for each group/clone to connect the shortest path
    mst_tree = mst(cluster_dist)
    sleep(0.5)

    edge_list = Tree()
    for c in tqdm(mst_tree, desc = 'Generating edge list '):
        G = nx.from_pandas_adjacency(mst_tree[c], create_using=nx.MultiDiGraph())
        G.edges(data=True)
        edge_list[c] = nx.to_pandas_edgelist(G)
    sleep(0.5)
    clone_ref = dict(out.metadata[clonekey])
    tmp_clone_tree = Tree()
    for x in out.metadata.index:
        tmp_clone_tree[clone_ref[x]][x].value = 1
    tmp_clone_tree2 = Tree()
    for x in tmp_clone_tree:
        tmp_clone_tree2[x] = list(tmp_clone_tree[x])

    tmp_clone_tree3 = Tree()
    for x in tmp_clone_tree2:
        tmp_ = pd.DataFrame(index = tmp_clone_tree2[x], columns = tmp_clone_tree2[x])
        tmp_ = pd.DataFrame(np.tril(tmp_) + 1, index = tmp_clone_tree2[x], columns = tmp_clone_tree2[x])
        tmp_.fillna(0, inplace = True)
        tmp_clone_tree3[x] = tmp_

    # here I'm using a temporary edge list to catch all cells that were identified as clones to forecfully link them up if they were clipped off during the mst step
    tmp_edge_list = Tree()
    for c in tqdm(tmp_clone_tree3, desc = 'Linking edges '):
        G = nx.from_pandas_adjacency(tmp_clone_tree3[c], create_using=nx.MultiDiGraph())
        G.edges(data=True)
        tmp_edge_list[c] = nx.to_pandas_edgelist(G)

    edge_listx = pd.concat([edge_list[x] for x in edge_list])
    edge_listx.index = [(s, t) for s, t in zip(edge_listx['source'],edge_listx['target'])]

    tmp_edge_listx = pd.concat([tmp_edge_list[x] for x in tmp_edge_list])
    tmp_edge_listx.index = [(s, t) for s, t in zip(tmp_edge_listx['source'], tmp_edge_listx['target'])]

    edge_list_final = edge_listx.combine_first(tmp_edge_listx)

    for idx in edge_list_final.index:
        edge_list_final.at[idx, 'weight'] = tmp_totaldist.loc[idx[0], idx[1]]
    # return the edge list
    edge_list_final.reset_index(drop = True, inplace = True)

    # and finally the vertex list which is super easy
    vertices = pd.DataFrame(out.metadata.index)

    # and now to actually generate the network
    edges = [tuple(e) for e in edge_list_final.values]
    edges_network = [tuple((e[0], e[1])) for e in edges]
    graph = igraph.Graph.Formula()
    graph.add_vertices(vertices['cell_id'])
    graph.add_edges(edges_network) # may take a very long time if there's too many

    if layout_option is not None:
        layout = graph.layout_fruchterman_reingold()
    else:
        layout = graph.layout(layout_option, *args, **kwds)
    for x in out.metadata.columns:
        graph.vs[x] = out.metadata[x]
    graph.es['width'] = [0.8/(int(e[2]) + 1) for e in edges]
    logg.info(' finished', time=start,
        deep=('Updated Dandelion object: \n'
        '   \'data\', contig-indexed clone table\n'
        '   \'metadata\', cell-indexed clone table\n'
        '   \'distance\', heavy and light chain distance matrices\n'
        '   \'edges\', network edges\n'
        '   \'layout\', network layout\n'
        '   \'graph\', network'))
    if self.__class__ == Dandelion:
        if self.germline is not None:
            germline_ = self.germline
        else:
            germline_ = None
        if self.threshold is not None:
            threshold_ = self.threshold
        else:
            threshold_ = None
        self.__init__(data = self.data, metadata = self.metadata, distance = dmat, edges = edge_list_final, layout = layout, graph = graph, germline = germline_, initialize = False)
        self.threshold = threshold_
    else:
        out = Dandelion(data = dat, distance = dmat, edges = edge_list_final, layout = layout, graph = graph, clone_key = clone_key)
        return(out)

def mst(mat):
    """
    Construct minimum spanning tree based on supplied matrix in dictionary.

    Parameters
    ----------
    mat : dict
        Dictionary containing numpy ndarrays.
    Returns
    ----------
        Dandelion `Tree` object holding DataFrames of constructed minimum spanning trees.
    """
    mst_tree = Tree()
    for c in mat:
        mst_tree[c] = pd.DataFrame(minimum_spanning_tree(np.triu(mat[c])).toarray().astype(int), index = mat[c].index, columns = mat[c].columns)
    return(mst_tree)

def transfer_network(self, dandelion, neighbors_key = None, rna_key = None, bcr_key = None):
    """
    Transfer data in `Dandelion` slots to `AnnData` object, updating the `.obs`, `.uns`, `.obsm` and `.raw` slots with metadata and network.

    Parameters
    ----------
    self : AnnData
        `AnnData` object
    dandelion : Dandelion
        `Dandelion` object
    neighbors_key : str, optional
        key for 'neighbors' slot in `.uns`
    rna_key : str, optional
        prefix for stashed RNA connectivities and distances.
    bcr_key : str, optional
        prefix for stashed BCR connectivities and distances.
    Returns
    ----------
        `AnnData` object with updated `.obs`, `.obsm` and '.obsp' slots with data from `Dandelion` object.

    """
    start = logg.info('Transferring network')
    if dandelion.edges is not None:
        G = nx.from_pandas_edgelist(dandelion.edges, create_using=nx.MultiDiGraph(), edge_attr='weight')
        distances = nx.to_pandas_adjacency(G, dtype = np.float32, weight='weight')
        connectivities = nx.to_pandas_adjacency(G, dtype = np.float32, weight=None)
        A = np.zeros(shape=(len(self.obs_names),len(self.obs_names)))
        df_connectivities = pd.DataFrame(A, index = self.obs_names, columns = self.obs_names)
        df_distances = pd.DataFrame(A, index = self.obs_names, columns = self.obs_names)
        print('converting matrices')
        df_connectivities.update(connectivities)
        df_distances.update(distances)

        df_connectivities_ = scipy.sparse.csr_matrix(df_connectivities.values, dtype = np.float32)
        df_distances_ = scipy.sparse.csr_matrix(df_distances.values, dtype = np.float32)

        print('Updating anndata slots')
        if neighbors_key is None:
            neighbors_key = "neighbors"
            rna_neighbors_key = 'rna_'+neighbors_key
            bcr_neighbors_key = 'rna_'+neighbors_key
            if rna_neighbors_key not in self.uns:
                self.uns[rna_neighbors_key] = self.uns[neighbors_key].copy()
            self.uns[bcr_neighbors_key] = {}
        if neighbors_key not in self.uns:
            raise ValueError("`edges=True` requires `pp.neighbors` to be run before.")

        if rna_key is None:
            r_connectivities_key = 'rna_connectivities'
            r_distances_key = 'rna_distances'
        else:
            r_connectivities_key = rna_key +'_connectivitites'
            r_distances_key = rna_key +'_distances'

        if bcr_key is None:
            b_connectivities_key = 'bcr_connectivities'
            b_distances_key = 'bcr_distances'
        else:
            b_connectivities_key = bcr_key +'_connectivitites'
            b_distances_key = bcr_key +'_distances'

        # stash_rna_connectivities:
        if r_connectivities_key not in self.obsp:
            self.obsp[r_connectivities_key] = self.obsp["connectivities"].copy()
            self.obsp[r_distances_key] = self.obsp["distances"].copy()

        # always overwrite the bcr slots
        self.obsp['connectivities'] = df_connectivities_.copy()
        self.obsp['distances'] = df_distances_.copy()
        self.obsp[b_connectivities_key] = self.obsp["connectivities"].copy()
        self.obsp[b_distances_key] = self.obsp["distances"].copy()

        self.uns[neighbors_key]['params'] = {'method':'bcr'}
        self.uns[bcr_neighbors_key] = self.uns[neighbors_key].copy()

    for x in dandelion.metadata.columns:
        self.obs[x] = pd.Series(dandelion.metadata[x])

    tmp = self.obs.copy()
    if dandelion.layout is not None:
        coord = pd.DataFrame(np.array(dandelion.layout), index = dandelion.metadata.index)
        for x in coord.columns:
            tmp[x] = coord[x]
        tmp[[1]] = tmp[[1]]*-1
        X_bcr = np.array(tmp[[0,1]], dtype = np.float32)
        self.obsm['X_bcr'] = X_bcr

    if (dandelion.edges is not None) and (dandelion.edges is not None):
        logg.info(' finished', time=start,
            deep=('updated `.obs` with `.metadata`\n'
                  'added to `.uns[\''+neighbors_key+'\']` and `.obsp`\n'
                  '   \'distances\', cluster-weighted adjacency matrix\n'
                  '   \'connectivities\', cluster-weighted adjacency matrix'))
    else:
        logg.info(' finished', time=start,
                deep=('updated `.obs` with `.metadata`\n'))

def define_clones(self, dist = None, action = 'set', model = 'ham', norm = 'len', doublets='drop', fileformat='airr', ncpu = None, dirs = None, outFilePrefix = None, key_added = None, verbose = False):
    """
    Find clones using changeo's `DefineClones.py <https://changeo.readthedocs.io/en/stable/tools/DefineClones.html>`__.

    Parameters
    ----------
    self : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones have been determined.
    dist : float, optional
        The distance threshold for clonal grouping. If None, the value will be retrieved from the Dandelion class .threshold slot.
    action : str
        Specifies how to handle multiple V(D)J assignments for initial grouping. Default is 'set'. The “first” action will use only the first gene listed. The “set” action will use all gene assignments and construct a larger gene grouping composed of any sequences sharing an assignment or linked to another sequence by a common assignment (similar to single-linkage).
    model : str
        Specifies which substitution model to use for calculating distance between sequences. Default is 'ham'. The “ham” model is nucleotide Hamming distance and “aa” is amino acid Hamming distance. The “hh_s1f” and “hh_s5f” models are human specific single nucleotide and 5-mer content models, respectively, from Yaari et al, 2013. The “mk_rs1nf” and “mk_rs5nf” models are mouse specific single nucleotide and 5-mer content models, respectively, from Cui et al, 2016. The “m1n_compat” and “hs1f_compat” models are deprecated models provided backwards compatibility with the “m1n” and “hs1f” models in Change-O v0.3.3 and SHazaM v0.1.4. Both 5-mer models should be considered experimental.
    norm : str
        Specifies how to normalize distances. Default is 'len'. 'none' (do not normalize), 'len' (normalize by length), or 'mut' (normalize by number of mutations between sequences).
    doublets : str
        Option to control behaviour when dealing with heavy chain 'doublets'. Default is 'drop'. 'drop' will filter out the doublets while 'count' will retain only the highest umi count contig.
    fileformat : str
        format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.
    ncpu : int, optional
        number of cpus for parallelization. Default is all available cpus.
    dirs : str, optional
        If specified, out file will be in this location.
    outFilePrefix : str, optional
        If specified, the out file name will have this prefix. None defaults to 'dandelion_define_clones'
    verbose : bool
        Whether or not to print the command used in terminal to call DefineClones.py. Default is False.
    Returns
    ----------
        `Dandelion` object with clone_id annotated in `.data` slot and `.metadata` initialized.
    """
    start = logg.info('Finding clones')
    if ncpu is None:
        nproc=multiprocessing.cpu_count()
    else:
        nproc=ncpu

    if key_added is None:
        clone_key = 'clone_id'
    else:
        clone_key = key_added

    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    else:
        dat = load_data(self)
    if os.path.isfile(str(self)):
        dat = load_data(self)
    dat_h = dat[dat['locus'] == 'IGH']
    dat_l = dat[dat['locus'].isin(['IGK', 'IGL'])]

    if os.path.isfile(str(self)):
        if dirs is None:
            tmpFolder = "{}/tmp".format(os.path.dirname(self))
            outFolder = "{}".format(os.path.dirname(self))
        else:
            tmpFolder = str(dirs).strip('/')+'/tmp'
            outFolder = str(dirs).strip('/')
    else:
        if dirs is None:
            tmpFolder = "dandelion/data/tmp"
            outFolder = "dandelion/data"
        else:
            tmpFolder = str(dirs).strip('/')+'/tmp'
            outFolder = str(dirs).strip('/')

    if not os.path.exists(tmpFolder):
        os.makedirs(tmpFolder)
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)

    if os.path.isfile(str(self)):
        h_file1 = "{}/{}_heavy-clone.tsv".format(tmpFolder, os.path.basename(self).split('.tsv')[0])
        h_file2 = "{}/{}_heavy-clone.tsv".format(outFolder, os.path.basename(self).split('.tsv')[0])
        l_file = "{}/{}_light.tsv".format(tmpFolder, os.path.basename(self).split('.tsv')[0])
        outfile = "{}/{}_clone.tsv".format(outFolder, os.path.basename(self).split('.tsv')[0])
    else:
        if outFilePrefix is not None:
            out_FilePrefix = outFilePrefix
        else:
            out_FilePrefix = 'dandelion_define_clones'
        h_file1 = "{}/{}_heavy-clone.tsv".format(tmpFolder, out_FilePrefix)
        h_file2 = "{}/{}_heavy-clone.tsv".format(outFolder, out_FilePrefix)
        l_file = "{}/{}_light.tsv".format(tmpFolder, out_FilePrefix)
        outfile = "{}/{}_clone.tsv".format(outFolder, out_FilePrefix)

    dat_h.to_csv(h_file1, sep = '\t', index = False)
    dat_l.to_csv(l_file, sep = '\t', index = False)

    if 'germline_alignment_d_mask' not in dat.columns:
        raise ValueError("Missing 'germline_alignment_d_mask' column in input file. Run create_germlines first.")

    if 'v_call_genotyped' in dat.columns:
        v_field = 'v_call_genotyped'
    else:
        v_field = 'v_call'

    if dist is None:
        if self.__class__ == Dandelion:
            if self.threshold is not None:
                dist_ = self.threshold
            else:
                raise ValueError('Threshold value in Dandelion object is None. Please run calculate_threshold first')
        else:
            raise ValueError('Distance value is None. Please provide a distance value (float)')
    else:
        dist_ = dist

    cmd = ['DefineClones.py',
            '-d', h_file1,
            '-o', h_file2,
            '--act', action,
            '--model', model,
            '--norm', norm,
            '--dist', str(dist_),
            '--nproc', str(nproc),
            '--vf', v_field]

    def clusterLinkage(cell_series, group_series):
        """
        Returns a dictionary of {cell_id : cluster_id} that identifies clusters of cells by analyzing their shared
        features (group_series) using single linkage.

        Arguments:
        cell_series (iter): iter of cell ids.
        group_series (iter): iter of group ids.

        Returns:
        dict:  dictionary of {cell_id : cluster_id}.
        """

        # assign initial clusters
        # initial_dict = {cluster1: [cell1], cluster2: [cell1]}
        initial_dict = {}
        for cell, group in zip(cell_series, group_series):
            try:
                initial_dict[group].append(cell)
            except KeyError:
                initial_dict[group] = [cell]

        # naive single linkage clustering (ON^2 best case, ON^3 worst case) ...ie for cells with multiple light chains
        # cluster_dict = {cluster1: [cell1, cell2]}, 2 cells belong in same group if they share 1 light chain
        while True:
            cluster_dict = {}
            for i, group in enumerate(initial_dict.keys()):
                cluster_dict[i] = initial_dict[group]
                for cluster in cluster_dict:
                    # if initial_dict[group] and cluster_dict[cluster] share common cells, add initial_dict[group] to cluster
                    if cluster != i and any(cell in initial_dict[group] for cell in cluster_dict[cluster]):
                        cluster_dict[cluster] = cluster_dict[cluster] + initial_dict[group]
                        del cluster_dict[i]
                        break
            # break if clusters stop changing, otherwise restart
            if len(cluster_dict.keys()) == len(initial_dict.keys()):
                break
            else:
                initial_dict = cluster_dict.copy()

        # invert cluster_dict for return
        assign_dict = {cell:k for k,v in cluster_dict.items() for cell in set(v)}

        return assign_dict

    def _lightCluster(heavy_file, light_file, out_file, doublets, fileformat):
        """
        Split heavy chain clones based on light chains

        Arguments:
        heavy_file (str): heavy chain input file.
        light_file (str): light chain input file.
        out_file (str): heavy chain output file.
        doublets (str): method for handling multiple heavy chains per cell. one of 'drop' or 'count'.
        format (str): file format. one of 'changeo' or 'airr'.
        """
        # Set column names
        if fileformat == 'changeo':
            cell_id = 'cell_id'
            clone_id = 'clone_id'
            v_call = 'v_call'
            j_call = 'j_call'
            junction_length = 'junction_length'
            umi_count = 'umicount'
        elif fileformat == 'airr':
            cell_id = 'cell_id'
            clone_id = 'clone_id'
            v_call = 'v_call'
            j_call = 'j_call'
            junction_length = 'junction_length'
            umi_count = 'umi_count'
        else:
            sys.exit("Invalid format %s" % fileformat)

        # read in heavy and light DFs
        heavy_df = pd.read_csv(heavy_file, dtype='object', na_values=['', 'None', 'NA'], sep='\t')
        light_df = pd.read_csv(light_file, dtype='object', na_values=['', 'None', 'NA'], sep='\t')

        # column checking
        expected_heavy_columns = [cell_id, clone_id, v_call, j_call, junction_length, umi_count]
        if set(expected_heavy_columns).issubset(heavy_df.columns) is False:
            raise ValueError("Missing one or more columns in heavy chain file: " + ", ".join(expected_heavy_columns))
        expected_light_columns = [cell_id, v_call, j_call, junction_length, umi_count]
        if set(expected_light_columns).issubset(light_df.columns) is False:
            raise ValueError("Missing one or more columns in light chain file: " + ", ".join(expected_light_columns))

        # Fix types
        heavy_df[junction_length] = heavy_df[junction_length].astype('int')
        light_df[junction_length] = light_df[junction_length].astype('int')

        # filter multiple heavy chains
        if doublets == 'drop':
            heavy_df = heavy_df.drop_duplicates(cell_id, keep=False)
            if heavy_df.empty is True:
                raise ValueError("Empty heavy chain data, after doublets drop. Are you combining experiments in a single file? If so, split your data into multiple files.")
        elif doublets == 'count':
            heavy_df[umi_count] = heavy_df[umi_count].astype('int')
            heavy_df = heavy_df.groupby(cell_id, sort=False).apply(lambda x: x.nlargest(1, umi_count))

        # transfer clone IDs from heavy chain df to light chain df
        clone_dict = {v[cell_id]:v[clone_id] for k, v in heavy_df[[clone_id, cell_id]].T.to_dict().items()}
        light_df = light_df.loc[light_df[cell_id].apply(lambda x: x in clone_dict.keys()), ]
        light_df[clone_id] = light_df.apply(lambda row: clone_dict[row[cell_id]], axis = 1)

        # generate a "cluster_dict" of CELL:CLONE dictionary from light df  (TODO: use receptor object V/J gene names)
        cluster_dict = clusterLinkage(light_df[cell_id],
                                    light_df.apply(lambda row:
                                                    getGene(row[v_call]) + ',' + \
                                                    getGene(row[j_call]) + ',' + \
                                                    str(row[junction_length]) + ',' + row[clone_id], axis=1))

        # add assignments to heavy_df
        heavy_df = heavy_df.loc[heavy_df[cell_id].apply(lambda x: x in cluster_dict.keys()), :]
        heavy_df[clone_id] = heavy_df[clone_id] + '_' + heavy_df.apply(lambda row: str(cluster_dict[row[cell_id]]), axis=1)

        # write heavy chains
        heavy_df.to_csv(out_file, sep='\t', index=False)
        return(heavy_df, light_df)

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd)

    h_df, l_df = _lightCluster(h_file2, l_file, outfile, doublets=doublets, fileformat=fileformat)

    h_df = load_data(h_df)
    # create a dictionary for cell_id : clone_id from h_df
    linked_clones = dict(zip(h_df['cell_id'], h_df['clone_id']))

    # create a clone_reference
    clone_ref = list(set(h_df['clone_id']))
    clone_ref = [c.split('_')[1] if c is not np.nan else c for c in clone_ref]
    l_df = load_data(l_df)

    for x in l_df.index:
        if l_df.loc[x, 'clone_id'] in clone_ref:
            l_df.at[x, 'clone_id'] = linked_clones[l_df.loc[x, 'cell_id']]
        else:
            try:
                l_df.at[x, 'clone_id'] = l_df.loc[x, 'cell_id']+'_notlinked'
            except:
                pass

    cloned_ = pd.concat([h_df, l_df])
    # transfer the new clone_id to the heavy + light file
    dat[str(clone_key)] = pd.Series(cloned_['clone_id'])

    if self.__class__ == Dandelion:
        if self.germline is not None:
            germline_ = self.germline
        else:
            germline_ = None
        if self.distance is not None:
            dist_ = self.distance
        else:
            dist_ = None
        if self.edges is not None:
            edge_ = self.edges
        else:
            edge_ = None
        if self.layout is not None:
            layout_ = self.layout
        else:
            layout_ = None
        if self.graph is not None:
            graph_ = self.graph
        else:
            graph_ = None
        if self.threshold is not None:
            threshold_ = self.threshold
        else:
            threshold_ = None

        if ('clone_id' in self.data.columns) and (clone_key is not None):
            self.__init__(data = dat, germline = germline_, distance = dist_, edges = edge_, layout = layout_, graph = graph_, initialize = True, retrieve = clone_key, split_heavy_light = False)
        elif ('clone_id' not in self.data.columns) and (clone_key is not None):
            self.__init__(data = dat, germline = germline_, distance = dist_, edges = edge_, layout = layout_, graph = graph_, initialize = True, clone_key = clone_key, retrieve = clone_key, split_heavy_light = False)
        else:
            self.__init__(data = dat, germline = germline_, distance = dist_, edges = edge_, layout = layout_, graph = graph_, initialize = True, clone_key = clone_key)
        self.threshold = threshold_
    else:
        if ('clone_id' in dat.columns) and (clone_key is not None):
            out = Dandelion(data = dat, retrieve = clonekey, split_heavy_light = False)
        elif ('clone_id' not in dat.columns) and (clone_key is not None):
            out = Dandelion(data = dat, clone_key = clone_key)
        else:
            out = Dandelion(data = dat)
        return(out)
    logg.info(' finished', time=start,
        deep=('Updated Dandelion object: \n'
        '   \'data\', contig-indexed clone table\n'
        '   \'metadata\', cell-indexed clone table\n'))

def quantify_clone_size(self, max_size = None, clone_key = None, key_added = None):
    """
    Quantifies size of clones

    Parameters
    ----------
    self : Dandelion
        `Dandelion` object
    max_size : int, optional
        The maximum size before value gets clipped. If None, the value will be returned as a numerical value.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata.
    key_added : str, optional
        Suffix to add to end of the output column.
    Returns
    ----------
        `Dandelion` object with clone size columns annotated in `.metadata` slot.
    """

    start = logg.info('Quantifying clone sizes')
    metadata_ = self.metadata.copy()

    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    clone_size = metadata_[str(clonekey)].value_counts()
    clone_group_size = metadata_[str(clonekey)+'_group'].value_counts()

    if max_size is not None:
        clone_size_ = clone_size.astype('object')
        clone_group_size_ = clone_group_size.astype('object')
        for i in clone_size.index:
            if clone_size.loc[i] >= max_size:
                clone_size_.at[i] = '>= '+ str(max_size)
        for i in clone_group_size.index:
            if clone_group_size.loc[i] >= max_size:
                clone_group_size_.at[i] = '>= '+ str(max_size)
        clone_size_ = clone_size_.astype('category')
        clone_group_size_ = clone_group_size_.astype('category')
    else:
        clone_size_ = clone_size.copy()
        clone_group_size_ = clone_group_size.copy()

    clone_size_dict = dict(clone_size_)
    clone_group_size_dict = dict(clone_group_size_)

    if key_added is None:
        self.metadata[str(clonekey)+'_size'] = pd.Series(dict(zip(metadata_.index, [clone_size_dict[c] for c in metadata_[str(clonekey)]])))
        self.metadata[str(clonekey)+'_group_size'] = pd.Series(dict(zip(metadata_.index, [clone_group_size_dict[c] for c in metadata_[str(clonekey)+'_group']])))
    else:
        self.metadata[str(clonekey)+'_size'+'_'+str(key_added)] = pd.Series(dict(zip(metadata_.index, [clone_size_dict[c] for c in metadata_[str(clonekey)]])))
        self.metadata[str(clonekey)+'_group_size'+'_'+str(key_added)] = pd.Series(dict(zip(metadata_.index, [clone_group_size_dict[c] for c in metadata_[str(clonekey)+'_group']])))
    logg.info(' finished', time=start,
        deep=('Updated Dandelion object: \n'
        '   \'metadata\', cell-indexed clone table'))
