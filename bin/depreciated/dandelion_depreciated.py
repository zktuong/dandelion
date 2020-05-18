#!/usr/bin/env python
import sys
if len(sys.argv) <= 2 :
    errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify one of the following options:\n\t1: Formatting\n\t2: Generate network\n\t2.1: Generate network (bcRep)\n\t2.2: Generate network (RBR)\n\t3: Plotting\n\t3.1: Plotting (bcRep)\n\t3.2: Plotting (RBR)\n\t3.2: Plotting (heavy chain only for option 2)\n'
    sys.exit(errormessage)
from time import time
from datetime import timedelta
dandelion_start_time = time()
import math
import os
import re
from collections import defaultdict, Iterable, Counter
import pandas as pd
import numpy as np
from distance import hamming
from Levenshtein import distance
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree
from subprocess import run
from Bio.Seq import translate
import multiprocessing
from joblib import Parallel, delayed
import seaborn as sns
from igraph import *

# from immcantation MakeDb.py
from presto.IO import readSeqFile, countSeqFile, printMessage, printError, printProgress, printLog
from Bio import SeqIO
from changeo.Gene import buildGermline
from changeo.IO import readGermlines, IgBLASTReader, ChangeoWriter, getOutputHandle, getFormatOperators

from changeo.Defaults import default_out_args
from changeo.Receptor import ChangeoSchema
from collections import OrderedDict
import csv

def checkEqual1(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return(True)
    return(all(first == rest for rest in iterator))

def Cluster_gini_index(clones_file, column = 'CLONE_IGHKL'):
    tab = pd.read_csv(clones_file, sep ='\t', dtype = 'object')
    clonesList = [int(x) for x in list(tab[column]) if str(x) != 'nan']
    freq = Counter(clonesList)
    points, vdf = VDF(list(freq.values()))
    gini = Get_Gini(points, vdf)
    return(gini)

def Cluster_i(fasta, tmp_file1, edge_length):
    cmd = "cd-hit -i "+fasta+" -o "+tmp_file1+" -c "+str(edge_length)+" -g 1 -d 80 -T 10 -M 0 -AL 40 -bak 1 -p 1"
    os.system(cmd)
    return()

def Cluster_sequences(seqs, cluster):
    clust_seq = Tree()
    for c in cluster:
        for id in cluster[c]:
            seqs1=[]
            if (id in seqs):
                s1=seqs[id]
                seqs1.append(s1)
            clust_seq[c][id] = seqs1
    return(clust_seq)

def extract(d, keys):
    return(dict((k, d[k]) for k in keys if k in d))

def Extract_edges(mst_tree, ld, L, cluster, file_out, file_out_h, threshold = 0):
    fh = open(file_out, 'w')
    fh.close()
    out = ''
    for c in mst_tree:
        nodes = []
        for i in cluster[c]:
            if type(i) == int:
                break
            nodes.append(i)
        A = mst_tree[c]
        if len(A) > 1:
            # convert to only lower diagonal
            A = A + A.T - np.diag(np.diag(A))
            A = np.tril(A)
            # keep a copy
            A2 = A
            # make edges unweighted
            A = np.where(A>0, 1, 0).astype(float)
            np.fill_diagonal(A, np.nan)
            B = ld[c]
            # add edge to non-diagonal identical sequences
            np.fill_diagonal(B, np.nan)
            A[np.tril(B==0)] = 1
            # also add edge to non-diagonal sequences that are different by 5% of the median length
            np.fill_diagonal(B, 0)
            A[np.tril(B<=(L[c] * threshold))] = 1
            np.fill_diagonal(A, 0)
            s, t = A.nonzero()
            ls = list(zip(s.tolist(), t.tolist()))
            for i in range(0, len(ls)):
                n1 = ls[i][0]
                n2 = ls[i][1]
                node1 = nodes[n1]
                node2 = nodes[n2]
                out = node1 + '\t' + node2 + '\t' + str(int(A2[ls[i]])) + '\n'
                Write_output(out, file_out)
                Write_output(out, file_out_h)
    return()

def Extract_edges_lightchain(file, clusters, barcode, contig, file_out, threshold = 0.01):
    all_seqs = Get_sequences(file)
    s_tree = Tree()
    hl_tree = Tree()
    for s in all_seqs:
        bc = s.split('_contig')[0]
        s_tree[bc][s].value = 1
    for bc1 in s_tree:
        if bc1 in barcode:
            hl_tree[bc1] = s_tree[bc1]
    hl_tree2 = Tree()
    for bc2 in hl_tree:
        for con in hl_tree[bc2]:
            if(con in contig):
                contig_list = list(hl_tree[bc2])
                hc_idx = contig_list.index(con)
                hc = contig_list[hc_idx]
                lc = contig_list[:hc_idx:]
                if len(lc) != 0:
                    hl_tree2[hc] = lc
    for w in hl_tree2:
        out = ''
        node1 = w
        for nn in hl_tree2[w]:
            node2 = ''.join(nn)
            out = str(node1) + '\t' + str(node2) + '\t' + str(200) +'\n'
            Write_output(out, file_out)

    for c in clusters:
        aa_list2 = []
        l2 = []
        light_chains = []
        for x in list(clusters[c]):
            light_chains.append(list(hl_tree2[x]))

        light_chains = list(flatten(light_chains))
        light_chains = list(set(light_chains))

        if len(light_chains) > 1:
            for nn2 in light_chains:
                aa_ = all_seqs[nn2]
                aa2 = translate(aa_, stop_symbol='@').split()
                aa_list2.append(aa2)
                l2.append(len(aa2[0]))
            lmean2 = math.floor(np.mean(l2))
            # prepare 2 dimensional array M x N (M entries (3) with N dimensions (1))
            tdarray2 = np.array(aa_list2).reshape(-1,1)
            d_mat2 = squareform(pdist(tdarray2,lambda x,y: distance(x[0],y[0])))
            A2 = d_mat2
            A2 = np.tril(A2)
            # keep a copy
            A3 = A2
            # make edges unweighted
            A2[np.tril(A2<=(lmean2 * threshold))] = 1
            A2[np.tril(A2>(lmean2 * threshold))] = 0
            np.fill_diagonal(A2, 0)
            source2, target2 = A2.nonzero()
            ls2 = list(zip(source2.tolist(), target2.tolist()))
            for i in range(0, len(ls2)):
                n1_2 = ls2[i][0]
                n2_2 = ls2[i][1]
                nn1_2 = light_chains[n1_2]
                nn2_2 = light_chains[n2_2]
                out_2 = str(nn1_2) + '\t' + str(nn2_2) + '\t' + str(int(A3[ls2[i]])) + '\n'
                Write_output(out_2, file_out)
    return()

def Extract_isotype(original_cr_anno, new_cr_anno):
    iso = {}
    iso1 = {}
    with open(original_cr_anno) as f:
        next(f)
        for line in f:
            line = line.strip()
            line = line.split(',')
            barcode = line[0]
            contig = line[2]
            c_gene = line[9]
            iso.setdefault(barcode, [])
            iso[barcode].append(c_gene)
            iso1.setdefault(barcode, [])
            iso1[barcode].append(contig)

        isotype = ['IGHA1', 'IGHA2', 'IGHD', 'IGHE', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHM']

        for barcode in iso:
            iso_match = which(pd.Series(iso[barcode]).isin(isotype).to_list())
            if len(iso_match) == 1:
                iso_ = []
                other = []
                for i in iso[barcode]:
                    if i in isotype:
                        iso_.append(i)
                    if i not in isotype:
                        other.append(i)
                if len(other) > 0:
                    iso[barcode] = [w.replace(other[0], iso_[0]) for w in iso[barcode]]

            if len(iso_match) > 1:
                x = [iso[barcode][i] for i in iso_match]
                if checkEqual1(x):
                    other = []
                    for i in iso[barcode]:
                        if i not in isotype:
                            other.append(i)
                    if len(other) > 0:
                        for o in range(0, len(other)):
                            iso[barcode] = [w.replace(other[o], x[0]) for w in iso[barcode]]
                else:
                    for i in iso[barcode]:
                        if i not in isotype:
                            iso[barcode][iso[barcode].index(i)] = 'ND'
            else:
                for i in iso[barcode]:
                    if i not in isotype:
                        iso[barcode][iso[barcode].index(i)] = 'ND'

    iso_list = [item for sublist in list(iso.values()) for item in sublist]
    contig_list = [item for sublist in list(iso1.values()) for item in sublist]

    cr_annotation = pd.read_csv(original_cr_anno)
    iso_df = pd.DataFrame(iso_list, contig_list)
    iso_df.rename(columns={0:'isotype'}, inplace=True)

    cr_annotation = cr_annotation.set_index('contig_id')
    new_anno = cr_annotation.join(iso_df)
    new_anno.to_csv(new_cr_anno)
    return()

def Extract_nodes(file_edges, seqs, file_out):
    fh2 = open(file_out, 'w')
    fh2.close()
    out = ''

    mapping = {0: 'IGHA1', 1: 'IGHA2', 2: 'IGHD', 3:'IGHE', 4:'IGHG1', 5:'IGHG2', 6:'IGHG3', 7:'IGHG4', 8:'IGHM'}
    fh1 = open(file_edges, 'r')
    v = []
    for l in fh1:
        l = l.strip()
        l = l.split()
        v.append(l)
    fh1.close()
    vertices = []
    for v1 in v:
        vertices.append(v1[0])
        vertices.append(v1[1])
    vertices = list(set(vertices))

    for i in vertices:
        iso = i.split('__')[1].split('_')
        iso_loc = np.array(iso).argmax()
        keys, inv = np.unique(iso_loc, return_inverse=True)
        vals = np.array([mapping[key] for key in keys])
        isotype = vals[inv]
        nodes = out + i + '\t' + str(1) + '\t' + ''.join(seqs[i]) + '\t' + ''.join(isotype) +'\n'
        Write_output(nodes, file_out)
    return()

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

def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el

def Filter_fasta(annotated_fasta, tab_file, filtered_fasta_out):
    seqs = Get_sequences(annotated_fasta)
    tab = pd.read_csv(tab_file, header=0, sep='\t')
    tab = tab[tab['FUNCTIONAL'] == 'T']
    filtered = extract(seqs, tab['SEQUENCE_ID'])
    fh = open(filtered_fasta_out, 'w')
    fh.close()
    for l in filtered:
        out = ''
        out=out+'>'+l+'\n'+filtered[l]+'\n'
        Write_output(out, filtered_fasta_out)
    return()

def filter_VDJ_fasta(tab_file, filtered_fasta_out):
    tab = pd.read_csv(tab_file, header=0, sep='\t')
    tab = tab[tab['FUNCTIONAL'] == 'T']
    tab.set_index('SEQUENCE_ID', inplace = True)
    seqs = dict(tab['SEQUENCE_VDJ'])
    for k, v in seqs.items():
        out = ''
        out = out + '>' + k + '\n' + v +'\n'
        Write_output(out, filtered_fasta_out)
    return()

def Filter_igh_fasta(annotated_fasta, tab_file, filtered_fasta_out):
    seqs = Get_sequences(annotated_fasta)
    tab_file = pd.read_csv(tab_file, header=0, sep='\t')
    tab_file = tab_file[(tab_file['FUNCTIONAL'] == 'T') & (tab_file['LOCUS'] == 'IGH')]
    filtered = extract(seqs, tab_file['SEQUENCE_ID'])
    fh = open(filtered_fasta_out, 'w')
    fh.close()
    for l in filtered:
        out = ''
        out=out+'>'+l+'\n'+filtered[l]+'\n'
        Write_output(out, filtered_fasta_out)
    return()

def Find_clones(tab_file, new_tab_file, clone_out, identity):
    ch = pd.read_csv(tab_file, sep = '\t', dtype = 'object')
    ch = ch.set_index('SEQUENCE_ID')
    ch1 = ch[(ch['FUNCTIONAL'].isin(['T'])) & (ch['LOCUS'].isin(['IGH']))]
    ch1 = ch1.dropna(subset=['JUNCTION_10X_AA'])

    # extracting Vgene, Jgene and junction length
    V = [re.sub('[*][0-9][0-9]', '', v) for v in ch1['V_CALL']]
    V = [', '.join(list(set(v.split(',')))) for v in V]
    J = [re.sub('[*][0-9][0-9]', '', j) for j in ch1['J_CALL']]
    J = [', '.join(list(set(j.split(',')))) for j in J]

    CDR3length = [len(str(l)) for l in ch1['JUNCTION_10X_AA']]
    junction = dict(zip(ch1.index, ch1['JUNCTION_10X_AA']))
    # combine into a dict and group IDs with same length, same V and J genes
    V_J_L = dict(zip(ch1.index, list(zip(V,J,CDR3length))))
    group = defaultdict(list)
    for key, val in sorted(V_J_L.items()):
        group[val].append(key)

    ids = Tree()
    seq = Tree()
    for g in group:
        id = group[g]
        j_aa = []
        idsl = []
        for i in id:
            j_aa.append(junction[i])
            idsl.append(i)
        seq[g] = j_aa
        ids[g] = idsl

    clones = Tree()
    for group in seq:
        unique = list(set(seq[group]))
        if len(unique) > 1:
            u_dict = dict(zip(range(0, len(unique)), unique))
        else:
            u_dict = dict({0:''.join(unique)})
        if len(unique) > 1:
            tdarray = np.array(unique).reshape(-1,1)
            d_mat = squareform(pdist(tdarray,lambda x,y: hamming(x[0],y[0])))
            tr = math.floor(group[2]*(1-identity))
            d_mat = np.tril(d_mat)
            # print(d_mat)
            np.fill_diagonal(d_mat, 0)
            indices_temp = [list(x) for x in np.tril_indices_from(d_mat)]
            indices = list(zip(indices_temp[0], indices_temp[1]))
            if len(indices) > 1:
                for pairs in indices:
                    if pairs[0] == pairs[1]:
                        indices.remove(pairs)
            indices_j = []
            for p in range(0, len(indices)):
                a1, b1 = indices[p]
                indices_j.append(u_dict[a1])
                indices_j.append(u_dict[b1])
            indices_j_f = list(set(indices_j))
            source, target = d_mat.nonzero()
            ls = list(zip(source.tolist(), target.tolist()))
            if len(source) == 0 & len(target) == 0:
                ls = list([(0,0)])
            dist = {}
            for l in ls:
                dist.update({l:d_mat[l]})
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
            # need a better way to try and extract the right pairs
            if len(dist) > 3:
                for d in dist:
                    if dist[d] < tr:
                        cm2.append(d)
                cm2 = list(set(cm2) ^ set(cm1))
            j_list1 = []
            if len(cm1) > 0:
                for l in range(0, len(cm1)):
                    a, b = cm1[l]
                    j_list1.append(u_dict[a])
                    j_list1.append(u_dict[b])
                j_list1 = list(set(j_list1))
            j_list2 = []
            if len(cm3) > 0:
                for l in range(0, len(cm2)):
                    a, b = cm2[l]
                    j_list2.append(u_dict[a])
                    j_list2.append(u_dict[b])
                j_list2 = list(set(j_list2))
            j_list3 = []
            if len(cm3) > 0:
                for l in range(0, len(cm3)):
                    a, b = cm3[l]
                    j_list3.append(u_dict[a])
                    j_list3.append(u_dict[b])
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
                clones[group+tuple(str(0))] = j_list1
            if len(j_list2) > 0:
                clones[group+tuple(str(1))] = j_list2
            if len(j_list3) > 0:
                for c in range(0, len(j_list3)):
                    clones[group+tuple(str(c+2))] = j_list3[c]
        else:
            clones[group+tuple(str(0))] = unique

    cid = Tree()
    for ck in clones.keys():
        ck2 = (ck[0], ck[1], ck[2])
        cid[ck][','.join(clones[ck])] = ids[ck2]

    clone_tree = Tree()
    # junction_tree = Tree()
    for x in cid:
        for y in cid[x]:
            keep = []
            # keep2 = []
            for z in cid[x][y]:
                id_junc = junction[z]
                if id_junc in y.split(','):
                    keep.append(z)
                    # keep2.append(id_junc)
            clone_tree[x][y] = keep
            # junction_tree[x][y] = keep2

    clone_df_dict = {}
    n = 1
    for t in dict(clone_tree):
        for k, v in dict(clone_tree[t]).items():
            clone_df_dict.update(dict({n:', '.join(v)}))
            n = n + 1

    ## assign clusters
    new_dict = {}
    for d in clone_df_dict:
        dd = clone_df_dict[d]
        ddd = dd.split(', ')
        for nd in ddd:
            new_dict.update(dict({nd:d}))

    ch_ = pd.DataFrame.from_dict(new_dict, orient = 'index', columns = ['CLONE_IGH'])
    ch_new = ch.merge(ch_.astype(str), 'outer', left_index = True, right_index = True)
    ch_new.reset_index(level=0, inplace=True)
    ch_new.rename(columns={'index':'SEQUENCE_ID'}, inplace=True)

    ch2 = ch_new[['SEQUENCE_ID', 'CLONE_IGH']]
    clones = ch2.set_index('SEQUENCE_ID').T.to_dict()
    hltree = Tree()
    for id in ch2['SEQUENCE_ID']:
        barcode = id.split('_contig')[0]
        hltree[barcode][id] = list(clones[id].values())

    for barcode in hltree:
        output = []
        for i in list(hltree[barcode].values()):
            output.append(i[0])
        if ~np.isnan(np.array(output).astype(float)).all():
            if np.isnan(np.array(output).astype(float)).any():
                clone = output[np.argwhere(~np.isnan(np.array(output).astype(float)))[0][0]]
                for id in hltree[barcode]:
                    hltree[barcode][id] = [clone]
            else:
                pass

    clone_ighkl = {}
    for i in hltree:
        clone_ighkl.update(dict(hltree[i]))

    clone_ighkl = pd.DataFrame.from_dict(clone_ighkl, orient ='index', columns = ['CLONE_IGHKL'])
    clone_ighkl.reset_index(level=0, inplace=True)
    clone_ighkl.rename(columns={'index':'SEQUENCE_ID'}, inplace=True)
    ch_new = ch_new.merge(clone_ighkl, 'outer', on='SEQUENCE_ID')
    ch_new.to_csv(new_tab_file, sep = '\t', index = False)

    ch_new_ss = ch_new[['SEQUENCE_ID', 'CLONE_IGH', 'CLONE_IGHKL']]
    barcode = []
    for id in ch_new_ss['SEQUENCE_ID']:
        barcode.append(id.split('_contig')[0])
    ch_out = pd.DataFrame(barcode, columns = ['BARCODE']).merge(ch_new_ss, left_index = True, right_index = True)
    ch_out.to_csv(clone_out, sep = '\t', header = None, index = None, na_rep = 'NA')

    return()

def Find_clones_bcRep(tab_file, rds_out, clone_out, new_tab_file):
    cmd = ['Rscript',
           '/Users/kt16/Documents/Clatworthy_scRNAseq/Ondrej/dandelion/bin/Find_clones_bcRep.R',
           '-i', tab_file,
           '-r', rds_out,
           '-o', clone_out,
           '-n', new_tab_file]
    print('Running command: %s\n' % (' '.join(cmd)))
    stdout = run(cmd)
    return(stdout)

def Format_annotation(id, iso_annotation, annotated_annotation):
    a = pd.read_csv(iso_annotation)
    iso_dict = {'IGHA1':'__1_0_0_0_0_0_0_0_0', 'IGHA2':'__0_1_0_0_0_0_0_0_0', 'IGHD':'__0_0_1_0_0_0_0_0_0', 'IGHE':'__0_0_0_1_0_0_0_0_0', 'IGHG1':'__0_0_0_0_1_0_0_0_0', 'IGHG2':'__0_0_0_0_0_1_0_0_0', 'IGHG3':'__0_0_0_0_0_0_1_0_0', 'IGHG4':'__0_0_0_0_0_0_0_1_0', 'IGHM':'__0_0_0_0_0_0_0_0_1', 'ND':'__0_0_0_0_0_0_0_0_0', 'ND':'__0_0_0_0_0_0_0_0_0', 'ND':'__0_0_0_0_0_0_0_0_0'}
    suffix = [iso_dict.get(iso) for iso in a['isotype']]
    new_contig_id = ''.join([str(a) + str(b) + '___' for a, b in list(zip(a['contig_id'],suffix))]).split('___')
    a['contig_id'] = [id+'_'+h for h in new_contig_id[:-1]]
    a.to_csv(annotated_annotation, index = False)
    return()

def Format_clones_RBR(cluster, tab_file, new_tab_file, clone_out):
    new_dict = {}
    for x in cluster:
        for y in cluster[x]:
            new_dict.update({y:x})
    for d in new_dict:
        new_dict[d] = int(new_dict[d])+1

    ch = pd.read_csv(tab_file, sep = '\t', dtype = 'object')
    if ('CLONE_RBR_IGH' in ch.columns) and ('CLONE_RBR_IGH' in ch.columns):
        ch.drop(columns = ['CLONE_RBR_IGH', 'CLONE_RBR_IGHKL'], inplace = True)
    ch = ch.set_index('SEQUENCE_ID')
    ch_ = pd.DataFrame.from_dict(new_dict, orient = 'index', columns = ['CLONE_RBR_IGH'])
    # ch_.CLONE_RBR_IGH = ch_.CLONE_RBR_IGH.astype(int)
    ch_new = ch.merge(ch_.astype(str), 'outer', left_index = True, right_index = True)
    ch_new.reset_index(level=0, inplace=True)
    ch_new.rename(columns={'index':'SEQUENCE_ID'}, inplace=True)

    ch2 = ch_new[['SEQUENCE_ID', 'CLONE_RBR_IGH']]
    clones = ch2.set_index('SEQUENCE_ID').T.to_dict()
    hltree = Tree()
    for id in ch2['SEQUENCE_ID']:
        barcode = id.split('_contig')[0]
        hltree[barcode][id] = list(clones[id].values())

    for barcode in hltree:
        output = []
        for i in list(hltree[barcode].values()):
            output.append(i[0])
        if ~np.isnan(np.array(output).astype(float)).all():
            if np.isnan(np.array(output).astype(float)).any():
                clone = output[np.argwhere(~np.isnan(np.array(output).astype(float)))[0][0]]
                for id in hltree[barcode]:
                    hltree[barcode][id] = [clone]
            else:
                pass

    clone_ighkl = {}
    for i in hltree:
        clone_ighkl.update(dict(hltree[i]))

    clone_ighkl_df = pd.DataFrame.from_dict(clone_ighkl, orient ='index', columns = ['CLONE_RBR_IGHKL'])
    clone_ighkl_df.reset_index(level=0, inplace=True)
    clone_ighkl_df.rename(columns={'index':'SEQUENCE_ID'}, inplace=True)
    ch_new = ch_new.merge(clone_ighkl_df, 'outer', on='SEQUENCE_ID')
    ch_new.to_csv(new_tab_file, sep = '\t', index = False)

    ch_new_ss = ch_new[['SEQUENCE_ID', 'CLONE_RBR_IGH', 'CLONE_RBR_IGHKL']]
    barcode = []
    for id in ch_new_ss['SEQUENCE_ID']:
        barcode.append(id.split('_contig')[0])
    ch_out = pd.DataFrame(barcode, columns = ['BARCODE']).merge(ch_new_ss, left_index = True, right_index = True)
    ch_out.to_csv(clone_out, sep = '\t', header = None, index = None, na_rep = 'NA')
    return()

def Format_header(fasta, annotation, annotated_fasta):
    a = pd.read_csv(annotation)
    iso_dict = {'IGHA1':'__1_0_0_0_0_0_0_0_0', 'IGHA2':'__0_1_0_0_0_0_0_0_0', 'IGHD':'__0_0_1_0_0_0_0_0_0', 'IGHE':'__0_0_0_1_0_0_0_0_0', 'IGHG1':'__0_0_0_0_1_0_0_0_0', 'IGHG2':'__0_0_0_0_0_1_0_0_0', 'IGHG3':'__0_0_0_0_0_0_1_0_0', 'IGHG4':'__0_0_0_0_0_0_0_1_0', 'IGHM':'__0_0_0_0_0_0_0_0_1', 'ND':'__0_0_0_0_0_0_0_0_0', 'ND':'__0_0_0_0_0_0_0_0_0', 'ND':'__0_0_0_0_0_0_0_0_0'}
    suffix = [iso_dict.get(iso) for iso in a['isotype']]
    new_header = ''.join([str(a) + str(b) + '___' for a, b in list(zip(a['contig_id'],suffix))]).split('___')
    a['suffix'] = new_id_header = [id+'_'+h for h in new_header[:-1]]
    a = a.set_index('contig_id')
    newheader_dict = a['suffix'].to_dict()
    fh = open(cellranger_fasta, 'r')
    seqs = {}
    for header, sequence in fasta_iterator(fh):
        seqs[newheader_dict[header]] =sequence
    fh.close()
    fh1 = open(annotated_fasta, 'w')
    fh1.close()
    out = ''
    for l in seqs:
        out = '>'+l+'\n'+seqs[l]+'\n'
        Write_output(out, annotated_fasta)
    return()

def Format_vertex_file(vertex_file, annotated_vertex_file, tab_file, meta_file, clone_column, pattern):
    vx_df = pd.read_csv(vertex_file, sep = '\t', header=None, na_values=[], keep_default_na = False, dtype = 'object')
    ctab = pd.read_csv(tab_file, sep = '\t', header=0, na_values=[], keep_default_na = False, dtype = 'object')
    meta = pd.read_csv(meta_file, sep = '\t', header=0, na_values=[], keep_default_na = False, dtype = 'object')
    donor = match('SANGERID', vx_df, meta, 0, 'SANGER SAMPLE ID', pattern)
    tissue = match('TISSUE', vx_df, meta, 0, 'SANGER SAMPLE ID', pattern)
    chain = match('LOCUS', vx_df, ctab, 0, 'SEQUENCE_ID', pattern = None)
    clone_igh = match(clone_column, vx_df, ctab, 0, 'SEQUENCE_ID', pattern = None)
    clone_ighkl = match(clone_column+'KL', vx_df, ctab, 0, 'SEQUENCE_ID', pattern = None)
    anno_vx_df = vx_df.assign(donor=donor, tissue = tissue, chain = chain, clone_igh = clone_igh, clone_ighkl = clone_ighkl)
    anno_vx_df.to_csv(annotated_vertex_file, sep = '\t', header = None, index = None, na_rep = 'NA')
    return()

def Generate_network_layout(vertex_file, edge_file, outdir, layout_file):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    nodes = pd.read_csv(vertex_file, sep = '\t', header=None, dtype = 'object')
    edge = pd.read_csv(edge_file, sep = '\t', header=None, dtype = 'object')
    edges = [tuple(e) for e in edge.values]
    edges = [e for e in edges if e[0] > e[1]]
    edges_list = [tuple((e[0], e[1])) for e in edges]
    g = Graph.Formula()
    g.add_vertices(nodes[0])
    g.add_edges(edges_list)
    g = g.simplify(combine_edges = 'first')
    layout = g.layout_graphopt(niter=800, node_charge = 0.00001, spring_constant = 3)
    out = pd.DataFrame(list(layout))
    out.to_csv(layout_file, sep = '\t', header = None, index = None)

    return(layout)

def Get_barcode_contig(cluster):
    bc = []
    contig = []
    for c in cluster:
        for id in cluster[c]:
            bc.append(id.split('_contig')[0])
            contig.append(id)
    bc = list(set(bc))
    return(bc, contig)

def Get_clusters(file):
    fh = open(file, "r")
    cluster = Tree()
    for l in fh:
        l = l.strip()
        l = l.split()
        id = l[1]
        cluster[l[2]][id].value = 1
    try:
        del cluster['NA']
    except:
        pass
    fh.close()
    return(cluster)

def Get_clusters_RBR(file):
    fh = open(file, "r")
    cluster = Tree()
    for l in fh:
        l = l.strip()
        l = l.split()
        id = l[2].replace("...", "")
        id = id.replace(">", "")
        cluster[l[0]][id].value = 1
    fh.close()
    return(cluster)

def Get_Gini(n,v):
    num_cores = multiprocessing.cpu_count()
    values=[]
    for i in range(0,len(n)):
        for j in range(0,v[i]):
            values.append(n[i])
    n = len(values)
    assert(n > 0), 'Empty list of values'
    sortedValues = sorted(values) #Sort smallest to largest
    cumm = [0]
    cumm.append(Parallel(n_jobs=num_cores)(delayed(sum_sorted_values)(i, sortedValues) for i in range(n)))
    cumm = list(flatten(cumm))
    LorenzPoints = [[], []]
    sumYs = 0           #Some of all y values
    robinHoodIdx = -1   #Robin Hood index max(x_i, y_i)
    for i in range(1, n + 2):
        x = 100.0 * (i - 1)/n
        y = 100.0 * (cumm[i - 1]/float(cumm[n]))
        LorenzPoints[0].append(x)
        LorenzPoints[1].append(y)
        sumYs += y
        maxX_Y = x - y
        if maxX_Y > robinHoodIdx: robinHoodIdx = maxX_Y
    giniIdx = 100 + (100 - 2 * sumYs)/n #Gini index
    return(giniIdx/100)

def Get_sequences(file, split_pattern = None):
    fh = open(file, "r")
    seqs = {}
    for header, sequence in fasta_iterator(fh):
        if split_pattern is not None:
            seqs[header.split(split_pattern)[0]] = sequence
        else:
            seqs[header] = sequence
    fh.close()
    return(seqs)

def Gini_indices(network, clones_file, id, out_file, column = 'CLONE_IGHKL'):
    fh = open(out_file, 'w')
    fh.close()
    out = ''
    vg = Vertex_knn_degree_gini_index(network)
    cg = Cluster_gini_index(clones_file, column)
    out = out + id + '\t' + str(vg) + '\t' + str(cg) + '\n'
    Write_output(out, out_file)
    return(vg, cg)

def igblastn(fasta, igblast_db, org='human'):
    nproc = multiprocessing.cpu_count()
    aux_db = os.path.join(igblast_db, 'optional_file', '%s_gl.aux' % org)
    v_db = os.path.join(igblast_db, 'database', 'imgt_'+org+'_ig_v')
    d_db = os.path.join(igblast_db, 'database', 'imgt_'+org+'_ig_d')
    j_db = os.path.join(igblast_db, 'database', 'imgt_'+org+'_ig_j')

    outsuffix = '_igblast.fmt7'
    out_file= fasta.split('.')[0]+outsuffix

    # igblastn command
    cmd = ['igblastn',
           '-query', os.path.abspath(fasta),
           '-out', os.path.abspath(out_file),
           '-num_threads', str(nproc),
           '-ig_seqtype', 'Ig',
           '-organism', org,
           '-auxiliary_data', str(aux_db),
           '-germline_db_V', str(v_db),
           '-germline_db_D', str(d_db),
           '-germline_db_J', str(j_db),
           '-outfmt', '7 std qseq sseq btop',
           '-domain_system', 'imgt']

    # set path to igblast database as environmental variable
    env = os.environ.copy()
    env['IGDATA'] = igblast_db
    print('Running command: %s\n' % (' '.join(cmd)))
    stdout = run(cmd, env=env)
    return(stdout)

def Levenshtein_dist_mat(clust_seq):
    ld = Tree()
    med_L = Tree()
    for c in clust_seq:
        aa_list = []
        l = []
        for i in clust_seq[c]:
            if type(i) == int:
                break
            s = clust_seq[c][i]
            # convert to amino acid
            aa = translate(''.join(s), stop_symbol='@').split()
            aa_list.append(aa)
            l.append(len(aa[0]))
        lmed = np.median(l)
        # prepare 2 dimensional array M x N (M entries (3) with N dimensions (1))
        tdarray = np.array(aa_list).reshape(-1,1)
        # calculate condensed distance matrix by wrapping the Levenshtein distance function
        d_mat = squareform(pdist(tdarray,lambda x,y: distance(x[0],y[0])))
        ld[c] = d_mat
        med_L[c] = lmed
    return(ld, med_L)

def match(return_column, query_df, target_df, query_column, target_column, pattern):
    if pattern is not None:
        if type(query_column) is int:
            if type(target_column) is int:
                x = [which(target_df.iloc[:, target_column] == j) for j in [re.split(pattern, i)[0] for i in list(query_df.iloc[:, query_column])]]
            if type(target_column) is str:
                x = [which(target_df.loc[:, target_column] == j) for j in [re.split(pattern, i)[0] for i in list(query_df.iloc[:, query_column])]]

        if type(query_column) is str:
            if type(target_column) is int:
                x = [which(target_df.iloc[:, target_column] == j) for j in [re.split(pattern, i)[0] for i in list(query_df.loc[:, query_column])]]
            if type(target_column) is str:
                x = [which(target_df.loc[:, target_column] == j) for j in [re.split(pattern, i)[0] for i in list(query_df.loc[:, query_column])]]
    else:
        if type(query_column) is int:
            if type(target_column) is int:
                x = [which(target_df.iloc[:, target_column] == j) for j in list(query_df.iloc[:, query_column])]
            if type(target_column) is str:
                x = [which(target_df.loc[:, target_column] == j) for j in list(query_df.iloc[:, query_column])]

        if type(query_column) is str:
            if type(target_column) is int:
                x = [which(target_df.iloc[:, target_column] == j) for j in list(query_df.loc[:, query_column])]
            if type(target_column) is str:
                x = [which(target_df.loc[:, target_column] == j) for j in list(query_df.loc[:, query_column])]

    result = []
    if type(return_column) is str:
        result = [''.join(target_df.loc[idx, return_column].replace(np.nan, 'NA', regex=True)) for idx in x]
    if type(return_column) is int:
        result = [''.join(target_df.iloc[idx, return_column].replace(np.nan, 'NA', regex=True)) for idx in x]
    return(result)

def mst(ld):
    mst_tree = Tree()
    for c in ld:
        mst_tree[c] = minimum_spanning_tree(ld[c]).toarray().astype(int)
    return(mst_tree)

def Plot_network_R(id, vertex_file, edge_file, layout_file, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    cmd = ['Rscript',
           '/Users/kt16/Documents/Clatworthy_scRNAseq/Ondrej/dandelion/bin/plotting_network.R',
           '-n', id,
           '-v', vertex_file,
           '-e', edge_file,
           '-l', layout_file,
           '-o', outdir]
    print('Running command: %s\n' % (' '.join(cmd)))
    stdout = run(cmd)
    return(stdout)

def Plot_network(vertex_file, edge_file, outdir, layout_file):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    nodes = pd.read_csv(vertex_file, sep = '\t', header=None, dtype = 'object')
    edge = pd.read_csv(edge_file, sep = '\t', header=None, dtype = 'object')

    if type(layout_file) != Layout:
        layout = pd.read_csv(layout_file, sep = '\t', header=None)
        layout = Layout(layout.values.tolist())
    else:
        layout = layout_file

    isotype_col_dict = {'IGHA1':'#4e79a7', 'IGHA2':'#f28e2b', 'IGHD':'#e15759', 'IGHE':'#76b7b2', 'IGHG1':'#59a14f', 'IGHG2':'#edc948', 'IGHG3':'#b07aa1', 'IGHG4':'#ff9da7', 'IGHM':'#9c755f', 'ND':'#c7c7c7'}
    tissue_col_dict = {'CAE':'#4E79A7', 'TCL':'#A0CBE8', 'SCL':'#F28E2B', 'TLN':'#FFBE7D', 'SPL':'#59A14F', 'LNG':'#8CD17D', 'MLN':'#B6992D', 'BMA':'#F1CE63', 'LIV':'#499894', 'DUO':'#86BCB6', 'ILE':'#E15759', 'KID':'#FF9D9A', 'OES':'#79706E', 'OME':'#BAB0AC', 'DKM':'#D37295', 'THY':'#FABFD2', 'TIL':'#B07AA1', 'SKM':'#D4A6C8'}
    donor_col_dict = {'A29':'#1f77b4', 'A31':'#ff7f0e', 'A35':'#2ca02c', 'A36':'#d62728', 'A37':'#9467bd'}
    chain_col_dict = {'IGH':'#c7c7c7', 'IGK':'#1F77B4', 'IGL':'#FF7F0E'}

    # file names
    outD = outdir + '/donor_network.pdf'
    outT = outdir + '/tissue_network.pdf'
    outC = outdir + '/cluster_network.pdf'
    outI = outdir + '/isotype_network.pdf'
    outL = outdir + '/locus_network.pdf'

    edges = [tuple(e) for e in edge.values]
    edges = [e for e in edges if e[0] > e[1]]
    edges_list = [tuple((e[0], e[1])) for e in edges]
    g = Graph.Formula()
    g.add_vertices(nodes[0])
    g.add_edges(edges_list)
    
    e_iso = []
    e_don = []
    e_tis = []
    e_chain = []
    e_cluster = []

    for e in edges:
        z = nodes[0].str.match(e[0])
        e_iso.append(nodes[3][z].to_string().split(' ')[4])
        e_don.append(nodes[4][z].to_string().split(' ')[4])
        e_tis.append(nodes[5][z].to_string().split(' ')[4])
        e_chain.append(nodes[6][z].to_string().split(' ')[4])
        e_cluster.append(nodes[8][z].to_string().split(' ')[4])
    e_anno = {0:e_iso, 1:e_don, 2:e_tis, 3:e_chain, 4:e_cluster}
    e_df = pd.DataFrame(e_anno)
    g.vs['isotype'] = nodes[3]
    g.vs['donor'] = nodes[4]
    g.vs['tissue'] = nodes[5]
    g.vs['chain'] = nodes[6]
    g.vs['cluster'] = nodes[8]
    g.es['isotype'] = e_df[0]
    g.es['donor'] = e_df[1]
    g.es['tissue'] = e_df[2]
    g.es['chain'] = e_df[3]
    g.es['cluster'] = e_df[4]
    g.es['width'] = [0.8/(int(e[2]) + 1) for e in edges]
    g = g.simplify(combine_edges = 'first')

    visual_style = {}
    visual_style['vertex_size'] = 1
    visual_style['vertex_frame_width'] = 0.1
    visual_style['vertex_label'] = g.vs['name']
    visual_style['vertex_label_size'] = 0
    # visual_style['edge_width'] = 0.05    
    visual_style['layout'] = layout
    visual_style["bbox"] = (300, 300)
    visual_style["margin"] = 20

    visual_style['edge_color'] = [isotype_col_dict[i] for i in g.es['isotype']]
    visual_style['vertex_color'] = [isotype_col_dict[i] for i in g.vs['isotype']]
    visual_style['vertex_frame_color'] = [isotype_col_dict[i] for i in g.vs['isotype']]
    plot(g, outI, **visual_style)

    visual_style['edge_color'] = [donor_col_dict[i] for i in g.es['donor']]
    visual_style['vertex_color'] = [donor_col_dict[i] for i in g.vs['donor']]
    visual_style['vertex_frame_color'] = [donor_col_dict[i] for i in g.vs['donor']]
    plot(g, outD, **visual_style)

    visual_style['edge_color'] = [tissue_col_dict[i] for i in g.es['tissue']]
    visual_style['vertex_color'] = [tissue_col_dict[i] for i in g.vs['tissue']]
    visual_style['vertex_frame_color'] = [tissue_col_dict[i] for i in g.vs['tissue']]
    plot(g, outT, **visual_style)

    visual_style['edge_color'] = [chain_col_dict[i] for i in g.es['chain']]
    visual_style['vertex_color'] = [chain_col_dict[i] for i in g.vs['chain']]
    visual_style['vertex_frame_color'] = [chain_col_dict[i] for i in g.vs['chain']]
    plot(g, outL, **visual_style)

    cluster_col_dict = dict(zip(list(set(nodes[8])), sns.color_palette("husl", len(list(set(nodes[8]))))))
    visual_style['edge_color'] = [cluster_col_dict[str(i)] for i in g.es['cluster']]
    visual_style['vertex_color'] = [cluster_col_dict[str(i)] for i in g.vs['cluster']]
    visual_style['vertex_frame_color'] = [cluster_col_dict[str(i)] for i in g.vs['cluster']]

    plot(g, outC, **visual_style)

    return(g)

def sum_sorted_values(i, sortedValues):
    res = sum(sortedValues[0:(i + 1)])
    return(res)

class Tree(defaultdict):
    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value

def VDF(n):
    points=sorted(n)
    vdf=[]
    for i in range(0,len(points)):
        vdf.append(n.count(points[i]))
    return(points, vdf)

def Vertex_knn_degree_gini_index(network):
    knn_degree = Graph.knn(network)[0]
    points, vdf = VDF(knn_degree)
    gini = Get_Gini(points, vdf)
    return(gini)

def which(self):
    try:
        self = list(iter(self))
    except TypeError as e:
        raise Exception("""'which' method can only be applied to iterables.
        {}""".format(str(e)))
    indices = [i for i, x in enumerate(self) if bool(x) == True]
    return(indices)

pd.Series.which = which

def Write_output(out, file):
    fh = open(file, "a")
    fh.write(out)
    fh.close()
    return()

### immcantation scripts for parsing igblast
cellranger_extended = ['cell', 'c_call', 'conscount', 'umicount',
                       'v_call_10x', 'd_call_10x', 'j_call_10x',
                       'junction_10x', 'junction_10x_aa', 'isotype']

def addGermline(receptor, references):

    __, germlines, __ = buildGermline(receptor, references)
    germline_seq = None if germlines is None else germlines['full']
    receptor.setField('germline_imgt', germline_seq)

    return(receptor)

def writeDb(records, fields, aligner_file, total_count, id_dict=None, annotations=None,
            partial=False, asis_id=True, writer=ChangeoWriter, out_file=None, out_args=default_out_args):

    def _open(x, f, writer=writer, out_file=out_file):
        if out_file is not None and x == 'pass':
            handle = open(out_file, 'w')
        else:
            handle = getOutputHandle(aligner_file,
                                     out_label='db-%s' % x,
                                     out_dir=out_args['out_dir'],
                                     out_name=out_args['out_name'],
                                     out_type=out_args['out_type'])
        return(handle, writer(handle, fields=f))

    # Function to convert fasta header annotations to changeo columns
    def _changeo(f, header):
        h = [ChangeoSchema.fromReceptor(x) for x in header if x.upper() not in f]
        f.extend(h)
        return(f)

    def _airr(f, header):
        h = [AIRRSchema.fromReceptor(x) for x in header if x.lower() not in f]
        f.extend(h)
        return(f)

    # Function to verify IMGT-gapped sequence and junction concur
    def _imgt_check(rec):
        try:  check = (rec.junction == rec.sequence_imgt[309:(309 + rec.junction_length)])
        except TypeError:  check = False
        return(check)

    # Function to check for valid records strictly
    def _strict(rec):
        valid = [rec.v_call and rec.v_call != 'None',
                 rec.j_call and rec.j_call != 'None',
                 rec.functional is not None,
                 rec.sequence_imgt,
                 rec.junction,
                 _imgt_check(rec)]
        return(all(valid))

    # Function to check for valid records loosely
    def _gentle(rec):
        valid = [rec.v_call and rec.v_call != 'None',
                 rec.d_call and rec.d_call != 'None',
                 rec.j_call and rec.j_call != 'None']
        return(any(valid))

    # Set writer class and annotation conversion function
    if writer == ChangeoWriter:
        _annotate = _changeo
    elif writer == AIRRWriter:
        _annotate = _airr
    else:
        printError('Invalid output writer.')

    # Set pass criteria
    _pass = _gentle if partial else _strict

    # Define log handle
    if out_args['log_file'] is None:
        log_handle = None
    else:
        log_handle = open(out_args['log_file'], 'w')

    # Initialize handles, writers and counters
    pass_handle, pass_writer = None, None
    fail_handle, fail_writer = None, None
    pass_count, fail_count = 0, 0
    start_time = time()

    # Validate and write outputf
    printProgress(0, total_count, 0.05, start_time=start_time)
    for i, record in enumerate(records, start=1):
        # Replace sequence description with full string, if required
        if id_dict is not None and record.sequence_id in id_dict:
            record.sequence_id = id_dict[record.sequence_id]

        # Parse sequence description into new columns
        if not asis_id:
            try:
                ann_raw = parseAnnotation(record.sequence_id)
                record.sequence_id = ann_raw.pop('ID')

                # Convert to Receptor fields
                ann_parsed = OrderedDict()
                for k, v in ann_raw.items():
                    ann_parsed[ChangeoSchema.toReceptor(k)] = v

                # Add annotations to Receptor and update field list
                record.setDict(ann_parsed, parse=True)
                if i == 1:  fields = _annotate(fields, ann_parsed.keys())
            except IndexError:
                # Could not parse pRESTO-style annotations so fall back to no parse
                asis_id = True
                printWarning('Sequence annotation format not recognized. Sequence headers will not be parsed.')

        # Add supplemental annotation fields
        # if _append_table is not None:
        #     record.setDict(_append_table(record.sequence_id), parse=True)
        if annotations is not None:
            record.setDict(annotations[record.sequence_id], parse=True)
            if i == 1:  fields = _annotate(fields, annotations[record.sequence_id].keys())

        # Count pass or fail and write to appropriate file
        if _pass(record):
            pass_count += 1
            # Write row to pass file
            try:
                pass_writer.writeReceptor(record)
            except AttributeError:
                # Open pass file and writer
                pass_handle, pass_writer = _open('pass', fields)
                pass_writer.writeReceptor(record)
        else:
            fail_count += 1
            # Write row to fail file if specified
            if out_args['failed']:
                try:
                    fail_writer.writeReceptor(record)
                except AttributeError:
                    # Open fail file and writer
                    fail_handle, fail_writer = _open('fail', fields)
                    fail_writer.writeReceptor(record)

        # Write log
        if log_handle is not None:
            log = OrderedDict([('ID', record.sequence_id),
                               ('V_CALL', record.v_call),
                               ('D_CALL', record.d_call),
                               ('J_CALL', record.j_call),
                               ('FUNCTIONAL', record.functional),
                               ('IMGT_PASS', _imgt_check(record))])
            printLog(log, log_handle)

        # Print progress
        printProgress(i, total_count, 0.05, start_time=start_time)

    # Print console log
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name) if pass_handle is not None else None
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    printLog(log)

    # Close file handles
    output = {'pass': None, 'fail': None}
    if pass_handle is not None:
        output['pass'] = pass_handle.name
        pass_handle.close()
    if fail_handle is not None:
        output['fail'] = fail_handle.name
        fail_handle.close()

    return(output)


def readCellRanger(cellranger_file, fields):
    # Mapping of 10X annotations to Receptor attributes
    cellranger_map = {'cell':  'barcode',
                      'c_call': 'c_gene',
                      'locus': 'chain',
                      'conscount': 'reads',
                      'umicount': 'umis',
                      'v_call_10x': 'v_gene',
                      'd_call_10x': 'd_gene',
                      'j_call_10x': 'j_gene',
                      'junction_10x': 'cdr3_nt',
                      'junction_10x_aa': 'cdr3',
                      'isotype': 'isotype'}

    # Function to parse individual fields
    def _parse(x):
        return('' if x == 'None' else x)

    # Generate annotation dictionary
    ann_dict = {}
    with open(cellranger_file) as csv_file:
        # Detect delimiters
        dialect = csv.Sniffer().sniff(csv_file.readline())
        csv_file.seek(0)
        # Read in annotation file
        csv_reader = csv.DictReader(csv_file, dialect=dialect)

        # Generate annotation dictionary
        for row in csv_reader:
            ann_dict[row['contig_id']] = {f: _parse(row[cellranger_map[f]]) for f in fields}

    return(ann_dict)

def getSeqDict(seq_file):
    seq_dict = SeqIO.to_dict(readSeqFile(seq_file), key_function=lambda x: x.description)

    return(seq_dict)

def parseIgBLAST(aligner_file, seq_file, repo, cellranger_file=None, partial=False, asis_id=True, asis_calls=False,
                 extended=False, format='changeo', out_file=None, out_args=default_out_args):
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MakeDB'
    log['ALIGNER'] = 'IgBLAST'
    log['ALIGNER_FILE'] = os.path.basename(aligner_file)
    log['SEQ_FILE'] = os.path.basename(seq_file)
    log['ASIS_ID'] = asis_id
    log['ASIS_CALLS'] = asis_calls
    log['PARTIAL'] = partial
    log['EXTENDED'] = extended
    printLog(log)

    start_time = time()
    printMessage('Loading files', start_time=start_time, width=20)

    # Count records in sequence file
    total_count = countSeqFile(seq_file)

    # Get input sequence dictionary
    seq_dict = getSeqDict(seq_file)

    # Create germline repo dictionary
    references = readGermlines(repo, asis=asis_calls)

    # Load supplementary annotation table
    if cellranger_file is not None:
        f = cellranger_extended if extended else cellranger_base
        annotations = readCellRanger(cellranger_file, fields=f)
    else:
        annotations = None

    printMessage('Done', start_time=start_time, end=True, width=20)

    # Check for IMGT-gaps in germlines
    if all('...' not in x for x in references.values()):
        printWarning('Germline reference sequences do not appear to contain IMGT-numbering spacers. Results may be incorrect.')

    # Define format operators
    try:
        __, writer, schema = getFormatOperators(format)
    except ValueError:
        printError('Invalid format %s.' % format)
    out_args['out_type'] = schema.out_type

    # Define output fields
    fields = list(schema.standard_fields)
    if extended:
        custom = IgBLASTReader.customFields(scores=True, regions=True, cdr3=False, schema=schema)
        fields.extend(custom)

    # Parse and write output
    with open(aligner_file, 'r') as f:
        parse_iter = IgBLASTReader(f, seq_dict, references, asis_calls=asis_calls)
        germ_iter = (addGermline(x, references) for x in parse_iter)
        output = writeDb(germ_iter, fields=fields, aligner_file=aligner_file, total_count=total_count,
                        annotations=annotations, partial=partial, asis_id=asis_id,
                        writer=writer, out_file=out_file, out_args=out_args)

    return(output)

### cmd options
id = sys.argv[1]
option = sys.argv[2]

outdir1 = id + '/dandelion/PLOT'
outdir2 = id + '/dandelion/DATA'
outdir3 = id + '/dandelion/CLONES'
outdir4 = id + '/dandelion/NETWORK'

if len(sys.argv) > 3 :
    plot_out = outdir1 + '/' + sys.argv[3]
meta_file = '/Users/kt16/Documents/Clatworthy_scRNAseq/Ondrej/dandelion/meta/PIP_sampleInfo_kt16.txt'
germline = ['/Users/kt16/Documents/Clatworthy_scRNAseq/Ondrej/dandelion/bcr_database/germlines/imgt/human/vdj/']
igdb = '/Users/kt16/Documents/Clatworthy_scRNAseq/Ondrej/dandelion/bcr_database/igblast'

seq_similarity = 0.95
junction_similarity = 0.85

cellranger_annotation = id + '/filtered_contig_annotations.csv'
cellranger_fasta = id + '/filtered_contig.fasta'

isotype_annotation = outdir2 + '/filtered_contig_annotations_isotype.csv'
annotated_annotation = outdir2 + '/filtered_contig_annotations_isotype_annotated.csv'
parsed_tab = outdir2 + '/filtered_contig_annotated_igblast_db-pass.tab'
clones_tab = outdir2 + '/filtered_contig_annotated_igblast_clone.tab'
clones_tab_RBR = outdir2 + '/filtered_contig_annotated_igblast_clone_RBR.tab'
clones_tab_bcrep = outdir2 + '/filtered_contig_annotated_igblast_clone_bcRep.tab'
igblast_out = outdir2 + '/filtered_contig_annotated_igblast.fmt7'
annotated_fasta = outdir2 + '/filtered_contig_annotated.fasta'
functional_fasta = outdir2 + '/filtered_contig_annotated_functional.fasta'
filtered_fasta_igh = outdir2 + '/filtered_contig_annotated_igh.fasta'
gini_out = outdir2 + '/gini_index.txt'
gini_out_bcrep = outdir2 + '/gini_index_bcRep.txt'
gini_out_RBR = outdir2 + '/gini_index_RBR.txt'

clones_out = outdir3 + '/filtered_contig_annotated_clones.txt'
clones_out_RBR = outdir3 + '/filtered_contig_annotated_RBR_clones.txt'
clones_out_bcrep = outdir3 + '/filtered_contig_annotated_bcRep_clones.txt'
clones_bcrep_rds = outdir3 + '/filtered_contig_annotated_bcRep_clones.RDS'
tmp_file = outdir3 + '/CDHIT_cluster_' + id + '.1'

file_vertices = outdir4 + '/Vertices_' + id + '.txt'
annotated_vertices = outdir4 + '/Annotated_vertices_' + id + '.txt'
file_edges = outdir4 + '/Edges_' + id + '.txt'
file_vertices_h = outdir4 + '/Vertices_heavy_' + id + '.txt'
annotated_vertices_h = outdir4 + '/Annotated_vertices_heavy_' + id + '.txt'
file_edges_h = outdir4 + '/Edges_heavy_' + id + '.txt'
file_vertices_bcrep = outdir4 + '/Vertices_' + id + '_bcRep.txt'
annotated_vertices_bcrep = outdir4 + '/Annotated_vertices_' + id + '_bcRep.txt'
file_edges_bcrep = outdir4 + '/Edges_' + id + '_bcRep.txt'
file_vertices_bcrep_h = outdir4 + '/Vertices_heavy_' + id + '_bcRep.txt'
annotated_vertices_bcrep_h = outdir4 + '/Annotated_vertices_heavy_' + id + '_bcRep.txt'
file_edges_bcrep_h = outdir4 + '/Edges_heavy_' + id + '_bcRep.txt'
file_vertices_RBR = outdir4 + '/Vertices_' + id + '_RBR.txt'
annotated_vertices_RBR = outdir4 + '/Annotated_vertices_' + id + '_RBR.txt'
file_edges_RBR = outdir4 + '/Edges_' + id + '_RBR.txt'
file_vertices_h_RBR = outdir4 + '/Vertices_heavy_' + id + '_RBR.txt'
annotated_vertices_h_RBR = outdir4 + '/Annotated_vertices_heavy_' + id + '_RBR.txt'
file_edges_h_RBR = outdir4 + '/Edges_heavy_' + id + '_RBR.txt'
layout_file = outdir4 + '/layout_'+ id + '.txt'
layout_file_bcrep = outdir4 + '/layout_'+ id + '_bcRep.txt'
layout_file_RBR = outdir4 + '/layout_'+ id + '_RBR.txt'

if(option == "1"):
    if not os.path.exists(outdir2):
        os.makedirs(outdir2)
    print('\nProcessing', '\t' + id + '\n')
    Extract_isotype(cellranger_annotation, isotype_annotation)
    Format_header(cellranger_fasta, isotype_annotation, annotated_fasta)
    Format_annotation(id, isotype_annotation, annotated_annotation)
    igblastn(annotated_fasta, igdb)
    parseIgBLAST(igblast_out, annotated_fasta, germline, annotated_annotation, extended=True)
    Filter_VDJ_fasta(parsed_tab, functional_fasta)
if(option == "2"):
    if not os.path.exists(outdir3):
        os.makedirs(outdir3)
    if not os.path.exists(outdir4):
        os.makedirs(outdir4)
    print('\nProcessing', '\t' + id + '\n')
    Find_clones(parsed_tab, clones_tab, clones_out, junction_similarity)
    cluster = Get_clusters(clones_out)
    seqs = Get_sequences(functional_fasta)
    clust_seqs = Cluster_sequences(seqs, cluster)
    ld, L = Levenshtein_dist_mat(clust_seqs)
    mst_tree = mst(ld)
    Extract_edges(mst_tree, ld, L, cluster, file_edges, file_edges_h)
    barcodes, contigs = Get_barcode_contig(cluster)
    Extract_edges_lightchain(functional_fasta, cluster, barcodes, contigs, file_edges)
    Extract_nodes(file_edges, seqs, file_vertices)
    Format_vertex_file(file_vertices, annotated_vertices, clones_tab, meta_file, 'CLONE_IGH', '_[ATCG][ATCG]')
if(option == "3"):
    if len(sys.argv) == 3 :
        errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
        sys.exit(errormessage)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)
    print('\nProcessing', '\t' + id + '\n')
    layout = Generate_network_layout(annotated_vertices, file_edges, plot_out, layout_file)
    Plot_network(annotated_vertices, file_edges, plot_out, layout)

if(option == "3.5"):
    if len(sys.argv) == 3 :
        errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
        sys.exit(errormessage)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)
    print('\nProcessing', '\t' + id + '\n')
    Plot_network(annotated_vertices, file_edges, plot_out, layout_file)

if(option == "3.5R"):
    if len(sys.argv) == 3 :
        errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
        sys.exit(errormessage)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)
    print('\nProcessing', '\t' + id + '\n')
    Plot_network_R(id, annotated_vertices, file_edges, layout_file, plot_out)

if(option == "4"):
    if len(sys.argv) == 3 :
        errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
        sys.exit(errormessage)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)
    print('\nProcessing', '\t' + id + '\n')
    network = Plot_network(annotated_vertices, file_edges, plot_out, layout_file)
    Gini_indices(network, clones_tab, id, gini_out, column = 'CLONE_IGHKL')

if(option == "2B"):
    if not os.path.exists(outdir3):
        os.makedirs(outdir3)
    if not os.path.exists(outdir4):
        os.makedirs(outdir4)
    print('\nProcessing', '\t' + id + '\n')
    try:
        Find_clones_bcRep(clones_tab, clones_bcrep_rds, clones_out_bcrep, clones_tab)
    except:
        Find_clones_bcRep(parsed_tab, clones_bcrep_rds, clones_out_bcrep, clones_tab_bcrep)
    cluster_bcrep = Get_clusters(clones_out_bcrep)
    seqs = Get_sequences(functional_fasta)
    clust_seqs = Cluster_sequences(seqs, cluster_bcrep)
    ld, L = Levenshtein_dist_mat(clust_seqs)
    mst_tree = mst(ld)
    Extract_edges(mst_tree, ld, L, cluster_bcrep, file_edges_bcrep, file_edges_bcrep_h)
    barcodes, contigs = Get_barcode_contig(cluster_bcrep)
    Extract_edges_lightchain(functional_fasta, cluster_bcrep, barcodes, contigs, file_edges_bcrep)
    Extract_nodes(file_edges_bcrep, seqs, file_vertices_bcrep)
    try:
        Format_vertex_file(file_vertices_bcrep, annotated_vertices_bcrep, clones_tab, meta_file, 'CLONE_bcRep_IGH','_[ATCG][ATCG]')
    except:
        Format_vertex_file(file_vertices_bcrep, annotated_vertices_bcrep, clones_tab_bcrep, meta_file, 'CLONE_bcRep_IGH','_[ATCG][ATCG]')

if(option == "3B"):
    if len(sys.argv) == 3 :
        errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
        sys.exit(errormessage)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)
    print('\nProcessing', '\t' + id + '\n')
    layout = Generate_network_layout(annotated_vertices_bcrep, file_edges_bcrep, plot_out, layout_file_bcrep)

if(option == "3BG"):
    if len(sys.argv) == 3 :
        errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
        sys.exit(errormessage)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)
    print('\nProcessing', '\t' + id + '\n')
    layout = Generate_network_layout(annotated_vertices_bcrep, file_edges_bcrep, plot_out, layout_file_bcrep)
    network = Plot_network(annotated_vertices_bcrep, file_edges_bcrep, plot_out, layout)

if(option == "3B_R"):
    if len(sys.argv) == 3 :
        errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
        sys.exit(errormessage)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)
    print('\nProcessing', '\t' + id + '\n')
    network = Plot_network(annotated_vertices_bcrep, file_edges_bcrep, plot_out, layout)
    Plot_network_R(id, annotated_vertices_bcrep, file_edges_bcrep, plot_out)

if(option == "4B"):
    if len(sys.argv) == 3 :
        errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
        sys.exit(errormessage)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)
    print('\nProcessing', '\t' + id + '\n')
    network = Plot_network(annotated_vertices_bcrep, file_edges_bcrep, plot_out, layout_file_bcrep)
    Gini_indices(network, clones_tab, id, gini_out_bcrep, column = 'CLONE_bcRep_IGHKL')

if(option == "2RBR"):
    if not os.path.exists(outdir3):
        os.makedirs(outdir3)
    if not os.path.exists(outdir4):
        os.makedirs(outdir4)
    print('\nProcessing', '\t' + id + '\n')
    Filter_igh_fasta(annotated_fasta, parsed_tab, filtered_fasta_igh)
    Cluster_i(filtered_fasta_igh, tmp_file, seq_similarity)
    cluster_RBR = Get_clusters_RBR(tmp_file+'.bak.clstr')
    try:
        Format_clones_RBR(cluster_RBR, clones_tab, clones_tab, clones_out_RBR)
    except:
        Format_clones_RBR(cluster_RBR, parsed_tab, clones_tab_RBR, clones_out_RBR)
    seqs_igh = Get_sequences(filtered_fasta_igh)
    clust_seqs_RBR = Cluster_sequences(seqs_igh, cluster_RBR)
    ld_igh, L_igh = Levenshtein_dist_mat(clust_seqs_RBR)
    mst_tree_full = mst(ld_igh)
    Extract_edges(mst_tree_full, ld_igh, L_igh, cluster_RBR, file_edges_RBR, file_edges_h_RBR)
    barcodes, contigs = Get_barcode_contig(cluster_RBR)
    Extract_edges_lightchain(functional_fasta, cluster_RBR, barcodes, contigs, file_edges_RBR)
    seqs = Get_sequences(functional_fasta)
    Extract_nodes(file_edges_RBR, seqs, file_vertices_RBR)
    try:
        Format_vertex_file(file_vertices_RBR, annotated_vertices_RBR, clones_tab, meta_file, 'CLONE_RBR_IGH', '_[ATCG][ATCG]')
    except:
        Format_vertex_file(file_vertices_RBR, annotated_vertices_RBR, clones_tab_RBR, meta_file, 'CLONE_RBR_IGH', '_[ATCG][ATCG]')
if(option == "3RBR"):
    if len(sys.argv) == 3 :
        errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
        sys.exit(errormessage)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)
    print('\nProcessing', '\t' + id + '\n')
    layout = Generate_network_layout(annotated_vertices_RBR, file_edges_RBR, plot_out, layout_file_RBR)
    network = Plot_network(annotated_vertices_RBR, file_edges_RBR, plot_out, layout)

if(option == "3RBR_R"):
    if len(sys.argv) == 3 :
        errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
        sys.exit(errormessage)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)
    print('\nProcessing', '\t' + id + '\n')
    Plot_network_R(id, annotated_vertices_RBR, file_edges_RBR, layout_file_RBR, plot_out_RBR)

if(option == "4RBR"):
    if len(sys.argv) == 3 :
        errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
        sys.exit(errormessage)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)
    print('\nProcessing', '\t' + id + '\n')
    network = Plot_network(annotated_vertices_RBR, file_edges_RBR, plot_out, layout_file_RBR)
    Gini_indices(network, clones_tab_RBR, id, gini_out_RBR, column = 'CLONE_RBR_IGHKL')

# if(option == "2.3"):
#     if not os.path.exists(outdir3):
#         os.makedirs(outdir3)
#     if not os.path.exists(outdir4):
#         os.makedirs(outdir4)
#     cluster = Get_clusters(clones_out)
#     seqs = Get_sequences(functional_fasta)
#     clust_seqs = Cluster_sequences(seqs, cluster)
#     ld, L = Levenshtein_dist_mat(clust_seqs)
#     mst_tree = mst(ld)
#     Extract_edges(mst_tree, ld, L, cluster, file_edges, file_edges_h)
#     Extract_nodes(file_edges, seqs, file_vertices)
#     Format_vertex_file(file_vertices, annotated_vertices, clones_tab, meta_file, '_[ATCG][ATCG]')
# if(option == "3.3"):
#     if len(sys.argv) == 3 :
#         errormessage = '\nUsage: dandelion.py <foldername> <option> <outname>\n\nPlease specify the out folder in outname\n'
#         sys.exit(errormessage)
#     if not os.path.exists(outdir1):
#         os.makedirs(outdir1)
#     seqs = Get_sequences(filtered_fasta_igh)
#     Extract_nodes(file_edges_h, seqs, file_vertices_h)
#     Format_vertex_file(file_vertices_h, annotated_vertices_h, clones_tab, meta_file, '_[ATCG][ATCG]')
#     network = Plot_network(annotated_vertices_h, file_edges_h, plot_out)
#     Gini_indices(network, clones_tab, id, gini_out_h)
#     Plot_network_R(id, annotated_vertices_h, file_edges_h, plot_out)

elapsed_time_secs = time() - dandelion_start_time
msg = "Execution took: %s secs (Wall clock time)\n" % timedelta(seconds=round(elapsed_time_secs))
print(msg)