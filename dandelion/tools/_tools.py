#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-05-13 23:22:18
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-05-22 15:55:39

import scanpy as sc
import pandas as pd
import numpy as np
from tqdm import tqdm
from ..utilities._misc import *
from collections import defaultdict
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
try:
    from scanpy import logging as logg
except ImportError:
    pass

def filter_bcr(data, adata, filter_bcr=True, filter_rna=True, filter_lightchains=True, filter_missing = True, outdir=None, outFilePrefix=None, filtered=False):
    dat = load_data(data)
    h = Tree()
    l = Tree()
    poor_qual = []
    h_doublet = []
    l_doublet = []

    locus_dict = dict(zip(dat['sequence_id'],dat['locus']))
    barcode = list(set(dat['cell_id']))

    bcr_check = Tree()
    for c in adata.obs_names:
        if c in barcode:
            bcr_check[c] = True
        else:
            bcr_check[c] = False
    adata.obs['has_bcr'] = pd.Series(dict(bcr_check))
    adata.obs['has_bcr'] = adata.obs['has_bcr'].astype('category')

    if 'v_call_genotyped' in dat.columns:
        v_dict = dict(zip(dat['sequence_id'], dat['v_call_genotyped']))
    else:
        v_dict = dict(zip(dat['sequence_id'], dat['v_call']))
    j_dict = dict(zip(dat['sequence_id'], dat['j_call']))

    for b in tqdm(barcode, desc = 'Marking barcodes with poor quality BCRs and BCR doublets'):
        hc_id = list(dat[(dat['cell_id'].isin([b])) & (dat['locus'] == 'IGH')]['sequence_id'])
        lc_id = list(dat[(dat['cell_id'].isin([b])) & (dat['locus'].isin(['IGK', 'IGL']))]['sequence_id'])
        h[b] = hc_id
        l[b] = lc_id
        # marking doublets defined by heavy chains
        if len(h[b]) > 1:
            h_doublet.append(b)
        # marking doublets defined by light chains
        if (len(h[b]) == 1) & (len(l[b]) > 1):
            l_doublet.append(b)
        # marking poor bcr quality, defined as those with only light chains, those
        # that were have conflicting assignment of locus and heavy/light V/J calls,
        # and also those that are missing either v or j calls
        if len(h[b]) < 1:
            poor_qual.append(b)
        if len(hc_id) > 0:
            v = v_dict[hc_id[0]]
            if 'IGH' not in v:
                poor_qual.append(b)
            j = j_dict[hc_id[0]]
            if 'IGH' not in j:
                poor_qual.append(b)
        if len(lc_id) > 0:
            v = v_dict[lc_id[0]]
            if 'IGH' in v:
                poor_qual.append(b)
            j = j_dict[lc_id[0]]
            if 'IGH' in j:
                poor_qual.append(b)

    poorqual = Tree()
    hdoublet = Tree()
    ldoublet = Tree()
    for c in tqdm(adata.obs_names, desc = 'Annotating in anndata obs slot '):
        if c in poor_qual:
            poorqual[c] = True
        else:
            poorqual[c] = False

        if c in h_doublet:
            hdoublet[c] = True
        else:
            hdoublet[c] = False

        if c in l_doublet:
            ldoublet[c] = True
        else:
            ldoublet[c] = False

    adata.obs['filter_bcr_quality'] = pd.Series(dict(poorqual))
    adata.obs['filter_bcr_quality'] = adata.obs['filter_bcr_quality'].astype('category')
    adata.obs['filter_bcr_heavy'] = pd.Series(dict(hdoublet))
    adata.obs['filter_bcr_heavy'] = adata.obs['filter_bcr_heavy'].astype('category')
    adata.obs['filter_bcr_light'] = pd.Series(dict(ldoublet))
    adata.obs['filter_bcr_light'] = adata.obs['filter_bcr_light'].astype('category')

    filter_ids = []
    if filter_bcr:
        if not filter_lightchains:
            filter_ids = list(set(h_doublet + poor_qual))
        else:
            filter_ids = list(set(h_doublet + l_doublet + poor_qual))

        if filter_rna:
            filter_ids = filter_ids + list(adata[adata.obs['filter_rna'] == True].obs_names)
            filter_ids = list(set(filter_ids))

        if filter_missing:
            for c in dat['cell_id']:
                if c not in adata.obs_names:
                    filter_ids.append(c)

        _dat = dat[~(dat['cell_id'].isin(filter_ids))]

        if os.path.isfile(str(data)):
            _dat.to_csv("{}/{}_filtered.tsv".format(os.path.dirname(data), os.path.basename(data).split('.tsv')[0]), sep = '\t', index = None)
        else:
            if filtered:
                outFile_prefix = 'filtered_contig'
            else:
                outFile_prefix = 'all_contig'

            if (outdir is None) & (outFilePrefix is not None):
                _dat.to_csv("{}/{}_filtered.tab".format('dandelion/data', str(outFilePrefix)), sep = '\t', index = None)
            elif (outdir is not None) & (outFilePrefix is None):
                _dat.to_csv("{}/{}_filtered.tab".format(str(outdir), outFile_prefix), sep = '\t', index = None)
            elif (outdir is None) & (outFilePrefix is None):
                _dat.to_csv("{}/{}_filtered.tab".format('dandelion/data', outFile_prefix), sep = '\t', index = None)
            elif (outdir is not None) & (outFilePrefix is None):
                _dat.to_csv("{}/{}_filtered.tab".format(str(outdir), outFile_prefix), sep = '\t', index = None)

        _adata = adata[~(adata.obs_names.isin(filter_ids))] # not saving the scanpy object because there's no need to at the moment

    return(_dat, _adata)

def find_clones(data, identity=0.85, outdir=None, clustering_by = None, by_alleles = None, outFilePrefix=None):
    dat = load_data(data)
    dat_heavy = dat[dat['locus'] == 'IGH']
    pd.set_option('mode.chained_assignment', None)

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
        raise ValueError("by_alleles only excepts boolean values or None, with None defaulting to False.")

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
        raise ValueError("clustering_by only excepts string values 'aa', 'nt' or None, with None defaulting to 'aa'.")
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
    for g in tqdm(junction_grp, desc = 'Finding clones '):
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
    dat_heavy['clone_id'] = pd.Series(clone_dict)
    hclone = dict(zip(dat_heavy['cell_id'], dat_heavy['clone_id']))
    hlclone = dict(zip(dat['sequence_id'], [hclone[c] for c in dat['cell_id']]))
    dat['clone_id'] = pd.Series(hlclone)
    # repeat this process for the light chains within each clone, but only for those with more than 1 light chains in a clone
    dat_light = dat[~(dat['locus'] == 'IGH')]
    # retrieve the J genes and J genes
    for c in tqdm(list(set(dat_light['clone_id'])), desc = 'Refining clone assignment based on light chain pairing '):
        dat_light_c = dat_light[dat_light['clone_id'] == c]
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
                raise ValueError("by_alleles only excepts boolean values or None, with None defaulting to False.")
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
                raise ValueError("clustering_by only excepts string values 'aa', 'nt' or None, with None defaulting to 'aa'.")

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
                dat.loc[renamed_clone_dict_light.keys(), 'clone_id'] = dat.loc[renamed_clone_dict_light.keys(), 'clone_id'] + '_' + pd.Series(renamed_clone_dict_light)

    if os.path.isfile(str(data)):
        dat.to_csv("{}/{}_clone.tsv".format(os.path.dirname(data), os.path.basename(data).split('.tsv')[0]), sep = '\t', index = False)
    else:
        if filtered:
            outFile_prefix = 'filtered_contig'
        else:
            outFile_prefix = 'all_contig'

        if (outdir is None) & (outFilePrefix is not None):
            dat.to_csv("{}/{}_clone.tsv".format('dandelion/data', str(outFilePrefix)), sep = '\t', index = None)
        elif (outdir is not None) & (outFilePrefix is None):
            dat.to_csv("{}/{}_clone.tsv".format(str(outdir), outFile_prefix), sep = '\t', index = None)
        elif (outdir is None) & (outFilePrefix is None):
            dat.to_csv("{}/{}_clone.tsv".format('dandelion/data', outFile_prefix), sep = '\t', index = None)
        elif (outdir is not None) & (outFilePrefix is None):
            dat.to_csv("{}/{}_clone.tsv".format(str(outdir), outFile_prefix), sep = '\t', index = None)

    return(dat)

def generate_network(data, distance_mode = None, clones_sep = None, layout_option = None, *args):
    """
    Extracting the necessary objects required for generating and plotting network.

    Arguments:
    data
        Dataframe in changeo/airr format after clones have been determined.
    clones_sep: tuple(int, str)
        A tuple containing how the clone groups should be extracted. None defaults to (0, '_')
    """
    start = logg.info('Generating network')
    dat = load_data(data)
    if 'clone_id' not in dat.columns:
        raise TypeError('Data does not contain clone information. Please run find_clones.')

    dat_h = dat[dat['locus'] == 'IGH']
    dat_l = dat[~(dat['locus'] == 'IGH')]
    if distance_mode is None or distance_mode is 'aa':
        seq_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h['sequence_alignment_aa'])))
        seq_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l['sequence_alignment_aa'])))
    elif distance_mode == 'nt':
        seq_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h['sequence_alignment'])))
        seq_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l['sequence_alignment'])))
    else:
        raise ValueError("distance_mode only excepts string values 'aa', 'nt' or None, with None defaulting to 'aa'.")

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
    total_dist = np.sum(dist_mat_list,axis=0)

    # Now, retrieve some 'meta data' for each cell that i can feed into obs lot in scanpy later
    # 1) clone
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
    # 3) clone_group
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
    metadata['clone_group'] = group
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

    # will return this object
    metadata = metadata[['clone_id', 'clone_group', 'isotype', 'lightchain', 'productive', 'heavychain_v', 'lightchain_v', 'heavychain_j', 'lightchain_j']]

    # generate edge list
    tmp_totaldist = pd.DataFrame(total_dist, index = metadata.index, columns = metadata.index)
    tmp_clusterdist = Tree()
    for i in metadata.index:
        cx = metadata.loc[i,'clone_group']
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
    clone_ref = dict(metadata['clone_id'])
    tmp_clone_tree = Tree()
    for x in metadata.index:
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
        edge_list_final.loc[idx, 'weight'] = tmp_totaldist.loc[idx[0], idx[1]]
    # return the edge list
    edge_list_final.reset_index(drop = True, inplace = True)

    # and finally the vertex list which is super easy
    vertices = pd.DataFrame(metadata.index)

    # and now to actually generate the network
    edges = [tuple(e) for e in edge_list_final.values]
    edges_network = [tuple((e[0], e[1])) for e in edges]
    graph = igraph.Graph.Formula()
    graph.add_vertices(vertices['cell_id'])
    graph.add_edges(edges_network) # may take a very long time if there's too many

    if layout_option is not None:
        layout = graph.layout_fruchterman_reingold()
    else:
        layout = graph.layout(layout_option, *args)
    for x in metadata.columns:
        graph.vs[x] = metadata[x]
    graph.es['width'] = [0.8/(int(e[2]) + 1) for e in edges]

    network = dandelion_network(data = dat, metadata = metadata, distance = d_mat, edges = edge_list_final, layout = layout, graph = graph)    
    logg.info(' finished', time=start,
        deep=('added to dandelion_network class object: \n'
        '   \'data\', contig-indexed clone table\n'
        '   \'metadata\', cell-indexed clone table\n'
        '   \'distance\', heavy and light chain distance matrices\n'
        '   \'edges\', network edges\n'
        '   \'layout\', network layout\n'
        '   \'graph\', network'))    
    return(network)

def mst(mat):
    mst_tree = Tree()
    for c in mat:
        mst_tree[c] = pd.DataFrame(minimum_spanning_tree(np.triu(mat[c])).toarray().astype(int), index = mat[c].index, columns = mat[c].columns)
    return(mst_tree)

def transfer_network(self, network, neighbors_key = None):
    start = logg.info('Transferring network')
    G = nx.from_pandas_edgelist(network.edges, create_using=nx.MultiDiGraph(), edge_attr='weight')
    distances = nx.to_pandas_adjacency(G, dtype = np.float32, weight='weight')
    connectivities = nx.to_pandas_adjacency(G, dtype = np.float32, weight=None)
    df_connectivities = pd.DataFrame(index = self.obs.index, columns = self.obs.index)
    df_distances = pd.DataFrame(index = self.obs.index, columns = self.obs.index)
    for x in connectivities.columns:
        df_connectivities[x] = pd.Series(connectivities[x])
    for x in distances.columns:
        df_distances[x] = pd.Series(distances[x])
    for x in df_distances.columns:
        df_distances[x] = df_distances[x].apply(lambda x: 5/(x + 1))
    df_connectivities.fillna(0, inplace = True)
    df_distances.fillna(0, inplace = True)
    df_connectivities_ = scipy.sparse.csr_matrix(df_connectivities.values, dtype = np.float32)
    df_distances_ = scipy.sparse.csr_matrix(df_distances.values, dtype = np.float32)
    
    if neighbors_key is None:
        neighbors_key = "neighbors"
    if neighbors_key not in self.uns:
        raise ValueError("`edges=True` requires `pp.neighbors` to be run before.")
    self.raw.uns = copy.deepcopy(self.uns)
    self.uns['neighbors']['connectivities'] = df_connectivities_
    self.uns['neighbors']['distances'] = df_distances_
    self.uns['neighbors']['params'] = {'method':'bcr'}

    for x in network.metadata.columns:
        self.obs[x] = pd.Series(network.metadata[x])

    tmp = self.obs.copy()
    coord = pd.DataFrame(np.array(network.layout), index = network.metadata.index)
    for x in coord.columns:
        tmp[x] = coord[x]
    tmp[[1]] = tmp[[1]]*-1
    X_bcr = np.array(tmp[[0,1]], dtype = np.float32)
    self.obsm['X_bcr'] = X_bcr
    logg.info(' finished', time=start,
        deep=('added to `.uns[\'neighbors\']`\n'
        '   \'distances\', cluster-weighted adjacency matrix\n'
        '   \'connectivities\', cluster-weighted adjacency matrix\n'
        'stored original .uns in .raw'))
    