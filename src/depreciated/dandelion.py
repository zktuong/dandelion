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
import functools
from functools import reduce

## this chunk here is not important yet
# def checkEqual1(iterator):
#     iterator = iter(iterator)
#     try:
#         first = next(iterator)
#     except StopIteration:
#         return(True)
#     return(all(first == rest for rest in iterator))

def Cluster_gini_index(clones_file, column = 'CLONE'):
    tab = pd.read_csv(clones_file, sep ='\t', dtype = 'object')
    clonesList = [int(x) for x in list(tab[column]) if str(x) != 'nan']
    freq = Counter(clonesList)
    points, vdf = VDF(list(freq.values()))
    gini = Get_Gini(points, vdf)
    return(gini)

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

def Gini_indices(network, clones_file, id, out_file, column = 'CLONE'):
    out = ''
    vg = Vertex_knn_degree_gini_index(network)
    cg = Cluster_gini_index(clones_file, column)
    out = out + id + '\t' + str(vg) + '\t' + str(cg) + '\n'
    Write_output(out, out_file)
    return(vg, cg)

def sum_sorted_values(i, sortedValues):
    res = sum(sortedValues[0:(i + 1)])
    return(res)

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

# the following code are essential
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

def Add_barcode_prefix(annotation, formatted_annotation, prefix = None):    
    anno = pd.read_csv(annotation, dtype = 'object')
    if prefix is not None:
        anno['barcode'] = [prefix+'_'+h.split('-')[0] for h in anno['barcode']]
        anno['contig_id'] = [prefix+'_'+h for h in anno['contig_id']]
        anno.index = [h.replace(prefix+'_', '') for h in anno['contig_id']]
    else:
        anno['barcode'] = [h for h in anno['barcode']]
        anno['contig_id'] = [h for h in anno['contig_id']]
        anno.index = [h for h in anno['contig_id']]
    
    anno.to_csv(formatted_annotation, index = None)
    return(anno)

def Format_fasta(fasta, annotation, out_fasta, prefix = None):
    if os.path.isfile(str(annotation)):
        anno = pd.read_csv(annotation, dtype = 'object', index_col = 0)
    else:
        anno = annotation
    if prefix is not None:
        anno.index = [h.replace(prefix+'_', '') for h in anno['contig_id']]
    else:
        anno.index = [h for h in anno['contig_id']]
    newheader_dict = anno['contig_id'].to_dict()
    fh = open(fasta, 'r')
    seqs = {}
    for header, sequence in fasta_iterator(fh):
        if header in newheader_dict:
            seqs[newheader_dict[header]] = sequence
    fh.close()
    fh1 = open(out_fasta, 'w')
    fh1.close()
    out = ''
    for l in seqs:
        out = '>'+l+'\n'+seqs[l]+'\n'
        Write_output(out, out_fasta)
        
def RunIgBlastn(fasta, igblast_db = None , org = 'human', loci = 'ig', format = 'blast', *args):
    env = os.environ.copy()
    if igblast_db is None:
        igblast_db = env['IGDATA']
    else:
        env['IGDATA'] = igblast_db
    
    cmd = ['AssignGenes.py', 'igblast',
           '-s', fasta,
           '-b', igblast_db,
           '--organism', org,
           '--loci', loci,
           '--format', format, *args]
    
    print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd, env=env)
    
def ParseIgBlast(igblast, fasta, germline, annotation, *args):    
    cmd = ['MakeDb.py', 'igblast',
            '-i', igblast,
            '-s', fasta,
            '-r', germline,
            '--10x', annotation,
            '--asis-id', *args]
    print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd)
    
def Transfer_isotype(changeo, col = 'C_CALL'):
    if os.path.isfile(str(changeo)):
        co = pd.read_csv(changeo, sep = '\t', dtype = 'object')
    else:
        co = changeo
    co['ISOTYPE'] = co[col]
    return(co)

class Tree(defaultdict):
    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value

def cleanNaNdict(d):
    return {
        k:v
        for k, v in d.items()
        if v is not np.nan
    }

def Check_barcodes(changeo, barcodes, regex_pattern = None, regex_split_selection = 0, barcode_columns = None):
    if os.path.isfile(str(changeo)):
        co = pd.read_csv(changeo, sep = '\t', dtype = 'object')
    else:
        co = changeo
    
    if type(barcodes) is dict:
        barcode_dict = barcodes
    elif (isinstance(barcodes, pd.DataFrame)) & (barcode_columns is not None):
        barcodes_ = barcodes
        if len(barcode_columns) == 2:
            barcode_dict = dict(zip(barcodes_[barcode_columns[0]], barcodes_[barcode_columns[1]]))
    elif (os.path.isfile(str(barcodes))) & (barcode_columns is not None):
        barcodes_ = pd.read_csv(barcodes, sep = '\t', dtype = 'object')
        if len(barcode_columns) == 2:
            barcode_dict = dict(zip(barcodes_[barcode_columns[0]], barcodes_[barcode_columns[1]]))

    barcode_dict = cleanNaNdict(barcode_dict)
    
    if (regex_pattern is not None) & (regex_split_selection is not None):
        keys = [re.split(regex_pattern, c)[regex_split_selection] for c in co['CELL']]
    
        for k in list(set(keys)):
            if k in barcode_dict.keys():
                co['CELL_'] = [c.replace(k, barcode_dict[k]) for c in co['CELL']]
    
    return(co, barcode_dict)

def Check_contigs(changeo, obs, obs_columns = None):
    if os.path.isfile(str(changeo)):
        co = pd.read_csv(changeo, sep = '\t', dtype = 'object')
    else:
        co = changeo
        
    if type(obs) is dict:
        filtered_dict = obs
    elif isinstance(obs, pd.DataFrame) & (obs_columns is not None):
        obs_ = obs.reset_index(drop = False)
        if len(obs_columns) == 2:
            filtered_dict = dict(zip(obs_[obs_columns[0]], obs_[obs_columns[1]]))
    elif (os.path.isfile(str(obs))) & (obs_columns is not None):
        obs_ = pd.read_csv(obs, sep = '\t', dtype = 'object')
        if len(obs_columns) == 2:
            filtered_dict = dict(zip(obs_[obs_columns[0]], obs_[obs_columns[1]]))
    
    newdict = {}
    nd = {}
    
    if 'CELL_' in co.columns:
        for c in co['CELL_']:
            try:
                matched = filtered_dict[c]
            except:
                matched = None
            newdict[c] = matched
        for k, v in co['CELL_'].items():
            nd[k] = newdict[v]
    else:
        for c in co['CELL']:
            try:
                matched = filtered_dict[c]
            except:
                matched = None
            newdict[c] = matched
        
        for k, v in co['CELL'].items():
            nd[k] = newdict[v]
    
    filterdf = pd.DataFrame.from_dict(nd, orient='index', columns = ['FILTER'])
    co['FILTER'] = filterdf
                
    return(co)

def Filter_changeo(changeo, changeo_filtered_pass, changeo_filtered_fail):
    if os.path.isfile(str(changeo)):
        co = pd.read_csv(changeo, sep = '\t', dtype = 'object')
    else:
        co = changeo    
    h = Tree()
    l = Tree()
    filter_ids = []
    
    co = co[(co['FUNCTIONAL'].isin(['T', 'TRUE', 'True', True])) & (co['IN_FRAME'].isin(['T', 'TRUE', 'True', True])) & (co['STOP'].isin(['F', 'FALSE', 'False', False]))]
    co = co[~(co['FILTER'].isin(['T', 'TRUE', 'True', True, np.nan]))]
    locus_dict = dict(zip(co['SEQUENCE_ID'],co['LOCUS']))
    barcode = list(set(co['CELL']))
    
    v_dict = dict(zip(co['SEQUENCE_ID'], co['V_CALL']))
    j_dict = dict(zip(co['SEQUENCE_ID'], co['J_CALL']))
    
    for b in barcode:
        hc_id = list(co[(co['CELL'].isin([b])) & (co['LOCUS'] == 'IGH')]['SEQUENCE_ID'])
        lc_id = list(co[(co['CELL'].isin([b])) & (co['LOCUS'].isin(['IGK', 'IGL']))]['SEQUENCE_ID'])
        h[b] = hc_id
        l[b] = lc_id
        if len(h[b]) != 1:
            filter_ids.append(b)
        if (len(h[b]) == 1) & (len(l[b]) > 1):            
            filter_ids.append(b)
        ## ok weird case. but sometimes the V/J call for the heavy locus is wrong?
        # should filter these away too
        if len(hc_id) > 0:
            v = v_dict[hc_id[0]]
            if 'IGH' not in v:
                filter_ids.append(b)
            j = j_dict[hc_id[0]]
            if 'IGH' not in j:
                filter_ids.append(b)        
                
    co1 = co[~(co['CELL'].isin(filter_ids))]
        
    co1.reset_index(drop = True, inplace = True)
    co1.to_csv(changeo_filtered_pass, sep = '\t', index = None)
    
    co2 = co[co['CELL'].isin(filter_ids)]
    h = Tree()
    l = Tree()
    filter_ids = []
    barcode = list(set(co2['CELL']))
    print('The following barcodes have more than 1 pair of Kappa AND/OR Lamda light chains:')    
    for b in barcode:
        hc_id = list(co2[(co2['CELL'].isin([b])) & (co2['LOCUS'] == 'IGH')]['SEQUENCE_ID'])
        lc_id = list(co2[(co2['CELL'].isin([b])) & (co2['LOCUS'].isin(['IGK', 'IGL']))]['SEQUENCE_ID'])
        h[b] = hc_id
        l[b] = lc_id
        if (len(h[b]) == 1) & (len(l[b]) > 1):
            filter_ids.append(b)
            print(b)
    co3 = co2[co2['CELL'].isin(filter_ids)]    
    co3.reset_index(drop = True, inplace = True)
    co3.to_csv(changeo_filtered_fail, sep = '\t', index = None)
    print('Refer to ambiguous file for more info.')    
    return(co1)

def Pad_seq(sequence):
    """
    Pad sequence to multiple of 3 with Ns
    """
    return {0: sequence, 1: sequence+'NN', 2: sequence+'N'}[len(sequence) % 3]

def cdhit(fasta, similarity = float, *args):
    cmd = ['cd-hit', 
           '-i', fasta, 
           '-o', fasta.replace('.fasta', '.1'), 
           '-c', str(similarity), 
           '-g', '1', 
           '-d', '80', 
           '-T', '10', 
           '-M', '0', 
           '-AL', '40', 
           '-bak', '1', 
           '-p', '1', *args]
    print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd)

def extract(d, keys):
    return(dict((k, d[k]) for k in keys if k in d))

    
def Split_fasta(fasta, changeo):
    if os.path.isfile(str(changeo)):
        co = pd.read_csv(changeo, sep = '\t', dtype = 'object')
    else:
        co = changeo    
    co_h = co[(co['FUNCTIONAL'].isin(['T', 'TRUE', 'True'])) & (co['IN_FRAME'].isin(['T', 'TRUE', 'True'])) & (co['STOP'].isin(['F', 'FALSE', 'False'])) & (co['LOCUS'] == 'IGH')]
    seqs = dict(zip(co_h['SEQUENCE_ID'], co_h['SEQUENCE_VDJ']))
    filtered_h = extract(seqs, co_h['SEQUENCE_ID'])
    fasta_h = fasta.replace('.fasta', '_vdj_heavy.fasta')
    fasta_l = fasta.replace('.fasta', '_vdj_light.fasta')

    fh = open(fasta_h, 'w')
    fh.close()
    for l in filtered_h:
        out = ''
        out=out+'>'+l+'\n'+filtered_h[l]+'\n'
        Write_output(out, fasta_h)

    co_l = co[(co['FUNCTIONAL'].isin(['T', 'TRUE', 'True'])) & (co['IN_FRAME'].isin(['T', 'TRUE', 'True'])) & (co['STOP'].isin(['F', 'FALSE', 'False'])) & (co['LOCUS'].isin(['IGK', 'IGL']))]
    seqs = dict(zip(co_l['SEQUENCE_ID'], co_l['SEQUENCE_VDJ']))
    filtered_l = extract(seqs, co_l['SEQUENCE_ID'])
    fh = open(fasta_l, 'w')
    fh.close()
    for l in filtered_l:
        out = ''
        out=out+'>'+l+'\n'+filtered_l[l]+'\n'
        Write_output(out, fasta_l)

def Retrieve_clusters(bakfile):
    fh = open(bakfile, "r")
    cluster = Tree()
    for l in fh:
        l = l.strip()
        l = l.split()
        id = l[2].replace("...", "")
        id = id.replace(">", "")
        cluster[l[0]][id].value = 1
    fh.close()
    return(cluster)

def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el

def Find_clones(changeo, changeoclone, identity, fasta = None, full_seq = False):
    start_time = time()
    if os.path.isfile(str(changeo)):
        co = pd.read_csv(changeo, sep = '\t', dtype = 'object')
    else:
        co = changeo
    co = co.set_index('SEQUENCE_ID')
    co1 = co[(co['FUNCTIONAL'].isin(['T', 'TRUE', 'True'])) & (co['IN_FRAME'].isin(['T', 'TRUE', 'True'])) & (co['STOP'].isin(['F', 'FALSE', 'False'])) & (co['LOCUS'] == 'IGH')]
    pd.set_option('mode.chained_assignment', None)
    co1['JUNCTION_AA'] = list(flatten([translate(Pad_seq(s), stop_symbol='@').split() for s in co1['JUNCTION']]))
        
    if (full_seq) & (fasta is not None):
        Split_fasta(fasta, changeo)
        cdhit(fasta.replace('.fasta', '_vdj_heavy.fasta'), identity)
        clones_full = Retrieve_clusters(fasta.replace('.fasta', '_vdj_heavy')+'.1.bak.clstr')
        new_dict = {}
        for x in clones_full:
            for y in clones_full[x]:
                new_dict.update({y:x})
        for d in new_dict:
            new_dict[d] = int(new_dict[d])+1        
    else:
        # extracting Vgene, Jgene and junction length
        V = [re.sub('[*][0-9][0-9]', '', v) for v in co1['V_CALL']]
        V = [', '.join(list(set(v.split(',')))) for v in V]
        J = [re.sub('[*][0-9][0-9]', '', j) for j in co1['J_CALL']]
        J = [', '.join(list(set(j.split(',')))) for j in J]

        CDR3length = [len(str(l)) for l in co1['JUNCTION_AA']]
        junction = dict(zip(co1.index, co1['JUNCTION_AA']))
        # combine into a dict and group IDs with same length, same V and J genes
        V_J_L = dict(zip(co1.index, list(zip(V,J,CDR3length))))
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
    
        for x in cid:
            for y in cid[x]:
                keep = []
                for z in cid[x][y]:
                    id_junc = junction[z]
                    if id_junc in y.split(','):
                        keep.append(z)
                clone_tree[x][y] = keep
    
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

    co_ = pd.DataFrame.from_dict(new_dict, orient = 'index', columns = ['CLONE'])
    co_new = co.merge(co_.astype(str), 'outer', left_index = True, right_index = True)  
    co_new.reset_index(level=0, inplace=True)
    co_new.rename(columns={'index':'SEQUENCE_ID'}, inplace=True)
    
    co2 = co_new[['SEQUENCE_ID', 'CLONE']]
    clones = co2.set_index('SEQUENCE_ID').T.to_dict()
    hltree = Tree()
    for id in co2['SEQUENCE_ID']:
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
            
    clone = {}
    for i in hltree:
        clone.update(dict(hltree[i]))
        
    clone = pd.DataFrame.from_dict(clone, orient ='index', columns = ['CLONE'])
    
    co_new.index = co_new['SEQUENCE_ID']
    co_new['CLONE'] = clone['CLONE']
    co_new = co_new[~(co_new['CLONE'].isin(['NA', 'NaN', 'nan', np.nan]))]
    co_new.to_csv(changeoclone, sep = '\t', index = False)
    
    elapsed_time_secs = time() - start_time
    msg = "Execution took: %s secs (Wall clock time)\n" % timedelta(seconds=round(elapsed_time_secs))
    print(msg)

    return(co_new)

def Find_clones_bcRep(changeo, changeoclone, identity):
    cmd = ['bcRep.R',
           '-c', changeo,
           '-o', changeoclone,
           '-i', str(identity)]
    run(cmd)
    
    co = pd.read_csv(changeoclone, sep = '\t', dtype = 'object')
    return(co) 

def immcantation_DefineClones(changeo_file, dist = 0.18, *args):
    start_time = time()
    cmd1 = ['ParseDb.py', 'select',
           '-d', changeo_file,
           '-f', 'LOCUS',
           '-u', 'IGH',
           '--logic', 'all',
           '--regex', '--outname',
           'TEST_heavy', *args]
    
    cmd2 = ['ParseDb.py', 'select',
           '-d', changeo_file,
           '-f', 'LOCUS',
           '-u', 'IG[LK]',
           '--logic', 'all',
           '--regex', '--outname',
           'TEST_light', *args]
    
    outPath = changeo_file.split('/')
    outPath3 = [re.sub('.*tab', 'TEST_heavy_parse-select.tab', h) for h in outPath]
    outPath3 = '/'.join(outPath3)
    cmd3 = ['DefineClones.py', 
            '-d', outPath3,
            '--act', 'set',
            '--model', 'ham',
            '--norm', 'len',
            '--dist', str(dist)]
    
    outPath4a = [re.sub('.*tab', 'TEST_heavy_parse-select_clone-pass.tab', h) for h in outPath]
    outPath4a = '/'.join(outPath4a)
    outPath4b = [re.sub('.*tab', 'TEST_light_parse-select.tab', h) for h in outPath]    
    outPath4b = '/'.join(outPath4b)
    outPath4c = [re.sub('.*tab', 'TEST_clone-pass.tab', h) for h in outPath]    
    outPath4c = '/'.join(outPath4c)
    
    cmd4 = ['light_cluster.py',
            '-d', outPath4a,
            '-e', outPath4b,
            '-o', outPath4c]
    
    print('Running command: %s\n' % (' '.join(cmd1)))
    run(cmd1)
    print('Running command: %s\n' % (' '.join(cmd2)))
    run(cmd2)
    print('Running command: %s\n' % (' '.join([str(c) for c in cmd3])))
    run(cmd3)
    print('Running command: %s\n' % (' '.join(cmd4)))
    run(cmd4)
    
    co = pd.read_csv(outPath4c, sep = '\t', dtype = 'object')
    elapsed_time_secs = time() - start_time
    msg = "Execution took: %s secs (Wall clock time)\n" % timedelta(seconds=round(elapsed_time_secs))
    print(msg)
    return(co)

def Get_clusters(changeoclone):
    if os.path.isfile(str(changeoclone)):
        co = pd.read_csv(changeoclone, sep = '\t', dtype = 'object')
    else:
        co = changeoclone
    co = co[(co['FUNCTIONAL'].isin(['T', 'TRUE', 'True'])) & (co['IN_FRAME'].isin(['T', 'TRUE', 'True'])) & (co['STOP'].isin(['F', 'FALSE', 'False'])) & (co['LOCUS'] == 'IGH')]
    clustdict = dict(zip(co['SEQUENCE_ID'], co['CLONE']))
    clust = Tree()
    for c in clustdict:
        clust[clustdict[c]][c].value = 1
    return(clust)

def Extract_sequences(changeoclone):
    if os.path.isfile(str(changeoclone)):
        co = pd.read_csv(changeoclone, sep='\t', index_col = 0, dtype = 'object')
    else:
        co = changeoclone
    co = co[(co['FUNCTIONAL'].isin(['T', 'TRUE', 'True'])) & (co['IN_FRAME'].isin(['T', 'TRUE', 'True'])) & (co['STOP'].isin(['F', 'FALSE', 'False']))]
    seqs = dict(zip(co['SEQUENCE_ID'], co['SEQUENCE_VDJ']))
    return(seqs)

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

def Levenshtein_dist_mat(clust_seq):
    ld = Tree()
    med_L = Tree()
    ld_ids = Tree()
    for c in clust_seq:
        aa_list = []
        l = []
        ids_list = []
        for i in clust_seq[c]:
            if type(i) == int:
                break
            # add some padding
            s = Pad_seq(''.join(clust_seq[c][i]))
            # convert gaps to N
            s = s.replace('-','N')
            # convert to amino acid
            aa = translate(s, stop_symbol='@').split()
            aa_list.append(aa)
            l.append(len(aa[0]))
            ids_list.append(i)
        lmed = np.median(l)
        # prepare 2 dimensional array M x N (M entries (3) with N dimensions (1)) 
        tdarray = np.array(aa_list).reshape(-1,1)
        # calculate condensed distance matrix by wrapping the Levenshtein distance function
        d_mat = squareform(pdist(tdarray,lambda x,y: distance(x[0],y[0])))
        ld[c] = d_mat
        med_L[c] = lmed
        ld_ids[c] = ids_list
    return(ld, med_L, ld_ids)

def mst(ld):
    mst_tree = Tree()
    for c in ld:
        mst_tree[c] = minimum_spanning_tree(ld[c]).toarray().astype(int)
    return(mst_tree)

def Extract_edges(mst_tree, ld, L, cluster, file_out, threshold = 0.01):
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
                
def Extract_nodes(file_edges, changeo, seqs, file_out):
    fh2 = open(file_out, 'w')
    fh2.close()
    out = ''
    
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
    
    if os.path.isfile(str(changeo)):
        co = pd.read_csv(changeo, sep = '\t', dtype = 'object')
    else:
        co = changeo
        
    co.index = co['SEQUENCE_ID']
    iso = co['ISOTYPE'].to_dict()
    
    for v in vertices:
        isotype = str(iso[v])
        nodes = out + v + '\t' + ''.join(seqs[v]) + '\t' + ''.join(isotype) +'\n'
        Write_output(nodes, file_out)

def Format_vertex_file(vertex_file, changeo_clone, meta_file, meta_column_key = None, meta_column_item = None, regex_pattern = None, regex_split_selection = 0, **kwargs):
    if os.path.isfile(str(changeo_clone)):
        co = pd.read_csv(changeo_clone, sep = '\t', header=0, na_values=[], keep_default_na = False, dtype = 'object')
    else:
        co = changeo_clone.reset_index(drop = True)
    vx_df = pd.read_csv(vertex_file, sep = '\t', header=None, na_values=[], keep_default_na = False, dtype = 'object')
    meta = pd.read_csv(meta_file, sep = '\t', header=0, na_values=[], keep_default_na = False, dtype = 'object')
    if (meta_column_key is not None) & (meta_column_item is not None):
        samp_dict = Construct_meta_dict(meta, meta_column_key, meta_column_item)
    co_dict = Construct_changeo_dict(co, **kwargs)
    
    tmp = Tree()
    if (meta_column_key is not None) & (meta_column_item is not None):
        for info in samp_dict:
            if (regex_pattern is not None) & (regex_split_selection is not None):
                tmp[info] = [samp_dict[info][re.split(regex_pattern, v)[regex_split_selection]] for v in vx_df[0]]
            else:
                tmp[info] = [samp_dict[info][v] for v in vx_df[0]]
    for info in co_dict:
        tmp[info] = [co_dict[info][v] for v in vx_df[0]]

    tmp_df = pd.DataFrame.from_dict(tmp)
    anno_vx_df = vx_df.join(tmp_df)
    anno_vx_df.to_csv(vertex_file.replace('.txt', '_annotated.txt'), sep = '\t', header = None, index = None, na_rep = 'NA')    
    
def Extract_dict(df, query):
    if isinstance(df, pd.DataFrame):
        query_ = df        
    elif os.path.isfile(str(df)):
        query_ = pd.read_csv(df, sep = '\t', dtype = 'object')
        
    if type(query) is tuple:
        if (type(query[0]) is list) & (type(query[1]) is not list):
            query_dict = dict(zip(query[0], query_[query[1]]))
        elif (type(query[1]) is list) & (type(query[0]) is not list):
            query_dict = dict(zip(query_[query[0]], query[1]))
        else:
            query_dict = dict(zip(query[0], query[1]))
    else:
        query_dict = dict(zip(query_[query[0]], query_[query[1]]))
        
    return(query_dict)

def Construct_meta_dict(meta, column_key, column_item):
    if os.path.isfile(str(meta)):
        info_ = pd.read_csv(meta, sep = '\t', dtype = 'object')
    else:
        info_ = meta
    dict_ = Tree()
    for i in column_item:
        dict_[i] = Extract_dict(info_, [column_key, i])
    return(dict_)

def Construct_changeo_dict(changeoclone, co_column_key = 'SEQUENCE_ID', co_column_item = ['LOCUS', 'CLONE', 'V_CALL', 'J_CALL', 'CELL'], strip_allele = None):
    if os.path.isfile(str(changeoclone)):
        co = pd.read_csv(changeoclone, sep = '\t', dtype = 'object')
    else:
        co = changeoclone.reset_index(drop = True)
    dict_ = Tree()
    if (strip_allele is None) | (not strip_allele):
        for i in co_column_item:
            dict_[i] = Extract_dict(co, (co_column_key, [re.sub('\*.[1234567890]', '', i_) for i_ in co[i]]))
    else:
        for i in co_column_item:
            dict_[i] = Extract_dict(co, [co_column_key, i])
    return(dict_)

def Format_combined_vertex_file(vertex_file, changeoclone, meta_file, meta_column_key = None, meta_column_item = None, barcode_dict = None, regex_pattern = None, regex_split_selection = 0, strip_allele = None, **kwargs):
    if os.path.isfile(str(changeoclone)):
        co = pd.read_csv(changeoclone, sep = '\t', header=0, na_values=[], keep_default_na = False, dtype = 'object')
    else:
        co = changeoclone.reset_index(drop = True)
    vx_df = pd.read_csv(vertex_file, sep = '\t', header=None, na_values=[], keep_default_na = False, dtype = 'object')
    meta = pd.read_csv(meta_file, sep = '\t', header=0, na_values=[], keep_default_na = False, dtype = 'object')
    if (meta_column_key is not None) & (meta_column_item is not None):
        samp_dict = Construct_meta_dict(meta, meta_column_key, meta_column_item)
    co_dict = Construct_changeo_dict(co, **kwargs)

    tmp = Tree()
    if (meta_column_key is not None) & (meta_column_item is not None):
        for info in samp_dict:
            if (regex_pattern is not None) & (regex_split_selection is not None):
                tmp[info] = [samp_dict[info][re.split(regex_pattern, v)[regex_split_selection]] for v in vx_df[0]]
            else:
                tmp[info] = [samp_dict[info][v] for v in vx_df[0]]
    for info in co_dict:
        tmp[info] = [co_dict[info][v] for v in vx_df[0]]

    tmp_df = pd.DataFrame.from_dict(tmp)
    anno_vx_df = vx_df.join(tmp_df)
    anno_vx_df['index'] = anno_vx_df['CELL']
    anno_vx_df.set_index('index', inplace = True)
    if barcode_dict is not None:
        if (regex_pattern is not None) & (regex_split_selection is not None):
            anno_vx_df['CELL_'] = [c.replace(re.split(regex_pattern, c)[regex_split_selection], barcode_dict[re.split(regex_pattern, c)[regex_split_selection]]) for c in anno_vx_df['CELL']]
        else:
            anno_vx_df['CELL_'] = [c.replace(c, barcode_dict[c]) for c in anno_vx_df['CELL']]

    cellbarcode = anno_vx_df['CELL']
    iso = Tree()
    vgene = Tree()
    jgene = Tree()
    isodict = Extract_dict(co, ['SEQUENCE_ID', 'ISOTYPE'])
    if (strip_allele is None) | (not strip_allele):
        vdict = Extract_dict(co, ('SEQUENCE_ID', [re.sub('\*.[1234567890]', '', i_) for i_ in co['V_CALL']]))
        jdict = Extract_dict(co, ('SEQUENCE_ID', [re.sub('\*.[1234567890]', '', i_) for i_ in co['J_CALL']]))
    else:
        vdict = Extract_dict(co, ['SEQUENCE_ID', 'V_CALL'])
        jdict = Extract_dict(co, ['SEQUENCE_ID', 'J_CALL'])
    for c in co['SEQUENCE_ID']:
        bc = c.split('-')[0]
        iso[bc][isodict[c]].value = 1
        vgene[bc][vdict[c]].value = 1
        jgene[bc][jdict[c]].value = 1
    ind2remove = ['IGKC', 'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7', np.nan]
    iso_assignment = Tree()
    v_assignment = Tree()
    j_assignment = Tree()
    for c in cellbarcode:
        iso_list = list(iso[c].keys())
        iso_list1 = [x for x in iso_list if x not in ind2remove]
        iso_assignment[c] = ''.join(iso_list1)
        v_list = list(vgene[c].keys())
        indices = [i for i, s in enumerate(v_list) if 'IGH' in s]
        if len(indices) > 0:
            v_list1 = v_list[indices[0]]
            v_list1 = list(set(v_list1.split(',')))
        else:
            v_list1 = ''
        if len(v_list1) > 1:
            v_assignment[c] = ','.join(v_list1)
        else:
            v_assignment[c] = ''.join(v_list1)
        j_list = list(jgene[c].keys())
        indices = [i for i, s in enumerate(j_list) if 'IGH' in s]
        if len(indices) > 0:
            j_list1 = j_list[indices[0]]
            j_list1 = list(set(j_list1.split(',')))
        else:
            j_list1 = ''
        if len(j_list1) > 1:
            j_assignment[c] = ','.join(j_list1)
        else:
            j_assignment[c] = ''.join(j_list1)
    iso_df = pd.DataFrame.from_dict(iso_assignment, orient='index', dtype = 'object', columns=['ISOTYPE'])
    v_df = pd.DataFrame.from_dict(v_assignment, orient='index', dtype = 'object', columns=['VGENE'])
    j_df = pd.DataFrame.from_dict(j_assignment, orient='index', dtype = 'object', columns=['JGENE'])
    
    tmp_dfs = [iso_df, v_df, j_df]
    tmp_df = reduce(lambda left,right: pd.merge(left,right,left_index = True, right_index = True), tmp_dfs)
    
    anno_vx_df = anno_vx_df.join(tmp_df)    
    anno_vx_df[2] = anno_vx_df['ISOTYPE']
    anno_vx_df['V_CALL'] = anno_vx_df['VGENE']
    anno_vx_df['J_CALL'] = anno_vx_df['JGENE']
    anno_vx_df.drop(['VGENE','JGENE'], axis = 1, inplace = True)
    anno_vx_df.to_csv(vertex_file.replace('.txt', '_annotated.txt'), sep = '\t', header = None, index = None, na_rep = 'NA')    
    anno_vx_df.reset_index(drop = False, inplace = True)
    if barcode_dict is not None:
        anno_vx_df['index'] = anno_vx_df['CELL_']
    anno_vx_df.reset_index(drop = True, inplace = True)
    anno_vx_df.set_index('index', inplace = True)
    anno_vx_df.drop([0,1,2], axis = 1, inplace = True)
    return(anno_vx_df)
    
def Generate_network_layout(vertex_file, edge_file, layout_file):
    nodes = pd.read_csv(vertex_file, sep = '\t', header=None, dtype = 'object')
    edge = pd.read_csv(edge_file, sep = '\t', header=None, dtype = 'object')
    edges = [tuple(e) for e in edge.values]
    edges = [e for e in edges if e[0] > e[1]]
    edges_list = [tuple((e[0], e[1])) for e in edges]
    g = Graph.Formula()
    g.add_vertices(nodes[0])
    g.add_edges(edges_list)
    # layout = g.layout_graphopt(niter=800, node_charge = 0.001, spring_constant = 3)
    layout = g.layout_fruchterman_reingold()
    out = pd.DataFrame(list(layout))
    out.to_csv(layout_file, sep = '\t', header = None, index = None)
    
    return(layout)

def Generate_network(vertex_file, edge_file):
    nodes = pd.read_csv(vertex_file, sep = '\t', header=None, dtype = 'object')
    edge = pd.read_csv(edge_file, sep = '\t', header=None, dtype = 'object')
     
    edges = [tuple(e) for e in edge.values]
    edges = [e for e in edges if e[0] > e[1]]
    edges_list = [tuple((e[0], e[1])) for e in edges]
    g = Graph.Formula()
    g.add_vertices(nodes[0])
    g.add_edges(edges_list)    
    
    iso_dict = dict(zip(nodes[0], nodes[2]))
    don_dict = dict(zip(nodes[0], nodes[3]))
    tis_dict = dict(zip(nodes[0], nodes[4]))
    locus_dict = dict(zip(nodes[0], nodes[5]))
    clone_dict = dict(zip(nodes[0], nodes[6]))

    e_iso = [iso_dict[e[0]] for e in edges]
    e_don = [don_dict[e[0]] for e in edges]
    e_tis = [tis_dict[e[0]] for e in edges]
    e_locus = [locus_dict[e[0]] for e in edges]
    e_clone = [clone_dict[e[0]] for e in edges]
    
    e_anno = {0:e_iso, 1:e_don, 2:e_tis, 3:e_locus, 4:e_clone}
    e_df = pd.DataFrame(e_anno)
    g.vs['isotype'] = nodes[2]
    g.vs['sample'] = nodes[3]
    g.vs['tissue'] = nodes[4]
    g.vs['locus'] = nodes[5]
    g.vs['clone'] = nodes[6]
    g.es['isotype'] = e_df[0]
    g.es['sample'] = e_df[1]
    g.es['tissue'] = e_df[2]
    g.es['locus'] = e_df[3]
    g.es['clone'] = e_df[4]
    g.es['width'] = [0.8/(int(e[2]) + 1) for e in edges]
    
    g = g.simplify(combine_edges = 'first')
    
    return(g)

def Plot_network(graph, layout, vertex_file, outdir, option):
    nodes = pd.read_csv(vertex_file, sep = '\t', header=None, dtype = 'object')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    if os.path.isfile(str(layout)):
        lyt = pd.read_csv(layout, sep = '\t', header=None, dtype = np.float32)
        lyt = Layout(lyt.values.tolist())
    else:
        lyt = layout
    
    isotype_col_dict = {'IGHA1':'#4e79a7', 'IGHA2':'#f28e2b', 'IGHD':'#e15759', 'IGHE':'#76b7b2', 'IGHG1':'#59a14f', 'IGHG2':'#edc948', 'IGHG3':'#b07aa1', 'IGHG4':'#ff9da7', 'IGHM':'#9c755f', 'NaN':'#f2f2f2', np.nan:'#f2f2f2', 'IGKC':'#f2f2f2', 'IGLC1':'#f2f2f2', 'IGLC2':'#f2f2f2', 'IGLC3':'#f2f2f2', 'IGLC4':'#f2f2f2', 'IGLC5':'#f2f2f2', 'IGLC6':'#f2f2f2', 'IGLC7':'#f2f2f2'}
    tissue_col_dict = {'CAE':'#4E79A7', 'TCL':'#A0CBE8', 'SCL':'#F28E2B', 'TLN':'#FFBE7D', 'SPL':'#59A14F', 'LNG':'#8CD17D', 'MLN':'#B6992D', 'BMA':'#F1CE63', 'LIV':'#499894', 'DUO':'#86BCB6', 'ILE':'#E15759', 'KID':'#FF9D9A', 'OES':'#79706E', 'OME':'#BAB0AC', 'DKM':'#D37295', 'THY':'#FABFD2', 'TIL':'#B07AA1', 'SKM':'#D4A6C8'}
    sample_col_dict = {'A29':'#1f77b4', 'A31':'#ff7f0e', 'A35':'#2ca02c', 'A36':'#d62728', 'A37':'#9467bd'}
    locus_col_dict = {'IGH':'#c7c7c7', 'IGK':'#1F77B4', 'IGL':'#FF7F0E'}
    
    g = graph
    
    # file names
    outD = outdir + '/sample_network.pdf'
    outT = outdir + '/tissue_network.pdf'
    outC = outdir + '/clone_network.pdf'
    outI = outdir + '/isotype_network.pdf'
    outL = outdir + '/locus_network.pdf'
    
    # set up visual style
    visual_style = {}
    visual_style['vertex_size'] = 1
    visual_style['vertex_frame_width'] = 0.1
    visual_style['vertex_label'] = g.vs['name']
    visual_style['vertex_label_size'] = 0
    # visual_style['edge_width'] = 0.05
    
    visual_style['layout'] = lyt
    visual_style['bbox'] = (300, 300)
    visual_style['margin'] = 20
    visual_style['inline'] = True
    
    if option == 'isotype':
        visual_style['edge_color'] = [isotype_col_dict[i] for i in g.es['isotype']]
        visual_style['vertex_color'] = [isotype_col_dict[i] for i in g.vs['isotype']]
        visual_style['vertex_frame_color'] = [isotype_col_dict[i] for i in g.vs['isotype']]
        p = plot(g, outI, **visual_style)
        
    if option == 'sample':
        visual_style['edge_color'] = [sample_col_dict[i] for i in g.es['sample']]
        visual_style['vertex_color'] = [sample_col_dict[i] for i in g.vs['sample']]
        visual_style['vertex_frame_color'] = [sample_col_dict[i] for i in g.vs['sample']]
        p = plot(g, outD, **visual_style)
    
    if option == 'tissue':
        visual_style['edge_color'] = [tissue_col_dict[i] for i in g.es['tissue']]
        visual_style['vertex_color'] = [tissue_col_dict[i] for i in g.vs['tissue']]
        visual_style['vertex_frame_color'] = [tissue_col_dict[i] for i in g.vs['tissue']]
        p = plot(g, outT, **visual_style)
    
    if option == 'locus':
        visual_style['edge_color'] = [locus_col_dict[i] for i in g.es['locus']]
        visual_style['vertex_color'] = [locus_col_dict[i] for i in g.vs['locus']]
        visual_style['vertex_frame_color'] = [locus_col_dict[i] for i in g.vs['locus']]
        p = plot(g, outL, **visual_style)
    
    if option == 'clone':
        clone_col_dict = dict(zip(list(set(nodes[6])), sns.color_palette("husl", len(list(set(nodes[6]))))))
        visual_style['edge_color'] = [clone_col_dict[str(i)] for i in g.es['clone']]
        visual_style['vertex_color'] = [clone_col_dict[str(i)] for i in g.vs['clone']]
        visual_style['vertex_frame_color'] = [clone_col_dict[str(i)] for i in g.vs['clone']]
        p = plot(g, outC, **visual_style)
    
    return(p)

def Plot_network_R(id, vertex_file, edge_file, layout_file, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    cmd = ['plotting_network.R',
            '-n', id,
            '-v', vertex_file,
            '-e', edge_file,
            '-l', layout_file,
            '-o', outdir]
    print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd)
    
def Get_barcode_contig(cluster):
    bc = []
    contig = []
    for c in cluster:
        for id in cluster[c]:
            bc.append(id.split('_contig')[0])
            contig.append(id)
    bc = list(set(bc))
    return(bc, contig)

def Extract_edges_lightchain(seqs, clusters, barcode, contig, file_out, threshold = 0):    
    all_seqs = seqs
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
                # add some padding
                aa_ = Pad_seq(aa_)
                # convert gaps to N
                aa_ = aa_.replace('-','N')
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
                
def Split_heavy_and_light(changeoclone):
    seq_out = Tree()
    if os.path.isfile(str(changeoclone)):
        co = pd.read_csv(changeoclone, sep='\t', index_col = 0, dtype = 'object')
    else:
        co = changeoclone
    co = co[(co['FUNCTIONAL'].isin(['T', 'TRUE', 'True'])) & (co['IN_FRAME'].isin(['T', 'TRUE', 'True'])) & (co['STOP'].isin(['F', 'FALSE', 'False']))]
    for clone in [x for x in list(set(co['CLONE'])) if str(x) != 'nan']:
        seq_out['heavy'][clone] = co[(co['CLONE'].isin([clone])) & (co['LOCUS'].isin(['IGH']))]['SEQUENCE_VDJ'].to_dict()
        seq_out['light'][clone] = co[(co['CLONE'].isin([clone])) & (co['LOCUS'].isin(['IGK', 'IGL']))]['SEQUENCE_VDJ'].to_dict()
    return(seq_out['heavy'], seq_out['light'])

def Calculate_combined_levenshstein_dist(changeoclone):
    seqh, seql = Split_heavy_and_light(changeoclone)
    ld_h, L_h, ld_ids_h = Levenshtein_dist_mat(seqh)
    ld_l, L_l, ld_ids_l = Levenshtein_dist_mat(seql)
    combined_ld = Tree()
    combined_ld_split = Tree()
    combined_ld_ids = Tree()
    for c in ld_ids_h:
        split_out = []
        HC = ld_h[c]
        LC = ld_l[c]
        HC_ID = [x.split('-')[0] for x in ld_ids_h[c]]
        LC_ID = [x.split('-')[0] for x in ld_ids_l[c]]
        if len(LC_ID) > 0:
            HC_df = pd.DataFrame(HC, index=HC_ID, columns=HC_ID)
            LC_df = pd.DataFrame(LC, index=LC_ID, columns=LC_ID)
            if len(HC_ID) == len(set(HC_ID)):
                if len(LC_ID) == len(set(LC_ID)):
                    if len(HC_ID) >= len(set(LC_ID)):
                        out = (LC_df.reindex_like(HC_df).fillna(0) + HC_df.fillna(0).fillna(0)).to_numpy()
                        split_out.append(out)
                        lids = ld_ids_l[c]
                        if len(HC_ID) > len(set(LC_ID)):
                            r = re.compile('|'.join(list(set(HC_ID) - set(LC_ID))))
                            ids = lids + list(filter(r.search, ld_ids_h[c]))
                        else:
                            ids = lids
                    elif len(HC_ID) < len(set(LC_ID)):
                        out = (HC_df.reindex_like(LC_df).fillna(0) + LC_df.fillna(0).fillna(0)).to_numpy()
                        split_out.append(out)
                        ids = ld_ids_l[c]
                else:                    
                    if len(HC_ID) >= len(set(LC_ID)):
                        indices = [i for i, s in enumerate([d == True for d in LC_df.index.duplicated(keep = False)]) if s]
                        for idx in indices:
                            LC_df1 = LC_df.copy()
                            LC_df1['idx'] = LC_df1.index
                            LC_df1.reset_index(drop=True, inplace=True)
                            LC_df1.drop(LC_df1.index[idx], inplace=True, axis=0)
                            LC_df1.drop(LC_df1.columns[idx], inplace=True, axis=1)
                            LC_df1.set_index('idx', inplace=True)
                            split_out.append((LC_df1.reindex_like(HC_df).fillna(0) + HC_df.fillna(0).fillna(0)).to_numpy())
                        out = sum(split_out)
                        lids = ld_ids_l[c]
                        if len(HC_ID) > len(set(LC_ID)):
                            r = re.compile('|'.join(list(set(HC_ID) - set(LC_ID))))
                            ids = lids + list(filter(r.search, ld_ids_h[c]))
                        else:
                            ids = lids
                    elif len(HC_ID) < len(set(LC_ID)):
                        indices = [i for i, s in enumerate([d == True for d in LC_df.index.duplicated(keep = False)]) if s]
                        for idx in indices:
                            LC_df1 = LC_df.copy()
                            LC_df1['idx'] = LC_df1.index
                            LC_df1.reset_index(drop=True, inplace=True)
                            LC_df1.drop(LC_df1.index[idx], inplace=True, axis=0)
                            LC_df1.drop(LC_df1.columns[idx], inplace=True, axis=1)
                            LC_df1.set_index('idx', inplace=True)
                            split_out.append((HC_df.reindex_like(LC_df1).fillna(0) + LC_df1.fillna(0).fillna(0)).to_numpy())
                        ids = ld_ids_l[c]
                        out = sum(split_out)
        else:
            out = HC
            split_out.append(out)
            ids = ld_ids_h[c]
              
        combined_ld[c] = out
        combined_ld_split[c] = split_out
        combined_ld_ids[c] = ids
    return(combined_ld, combined_ld_split, combined_ld_ids)

def Extract_combined_edges(comb_mst_tree, combined_ld, combined_ld_ids, file_out):
    fh = open(file_out, 'w')
    fh.close()
    out = ''
    for c in comb_mst_tree:
        nodes = []
        for i in combined_ld_ids[c]:
            if type(i) == int:
                break
            nodes.append(i)
        A = comb_mst_tree[c]
        if len(A) > 1:
            # convert to only lower diagonal
            A = A + A.T - np.diag(np.diag(A))
            A = np.tril(A)
            # keep a copy
            A2 = A
            # make edges unweighted
            A = np.where(A>0, 1, 0).astype(float)
            np.fill_diagonal(A, np.nan)
            B = combined_ld[c]
            # add edge to non-diagonal identical sequences
            np.fill_diagonal(B, np.nan)
            A[np.tril(B==0)] = 1
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

### cmd options
id = sys.argv[1]
option = sys.argv[2]

outdir1 = id + '/dandelion/data'
outdir2 = id + '/dandelion/figures'
outdir3 = id + '/dandelion/network'

if len(sys.argv) > 3 :
    plot_out = outdir2 + '/' + sys.argv[3]

if not os.path.exists(outdir1):
    os.makedirs(outdir1)
if not os.path.exists(outdir2):
    os.makedirs(outdir2)
if not os.path.exists(outdir3):
    os.makedirs(outdir3)

sampleInfo_file = '/Users/kt16/Documents/Clatworthy_scRNAseq/Ondrej/dandelion_files/meta/PIP_sampleInfo_kt16.txt'
obs_file = '/Users/kt16/Documents/Clatworthy_scRNAseq/Ondrej/dandelion_files/meta/pip_scanpy_obs.txt'

germline = '/Users/kt16/Documents/Github/dandelion/database/germlines/imgt/human/vdj/'
igdb = '/Users/kt16/Documents/Github/dandelion/database/igblast'

cellranger_annotation = id + '/filtered_contig_annotations.csv'
cellranger_fasta = id + '/filtered_contig.fasta'
filtered_fasta = outdir1 + '/filtered_contig_annotated_filtered.fasta'

annotation_ = outdir1 + '/filtered_contig_annotations_.csv'
filtered_annotation = outdir1 + '/filtered_contig_annotations_filtered.csv'
annotated_annotation = outdir1 + '/filtered_contig_annotations_filtered_annotated.csv'

parsed_tab = outdir1 + '/filtered_contig_annotated_igblast_db-pass.tab'
filtered_tab = outdir1 + '/filtered_contig_annotated_igblast_db_filtered.tab'
ambiguous_tab = outdir1 + '/filtered_contig_annotated_igblast_db_filtered-ambiguous.tab'
clones_tab = outdir1 + '/filtered_contig_annotated_igblast_db_filtered_clone.tab'

igblast_out = outdir1 + '/filtered_contig_annotated_igblast.fmt7'
annotated_fasta = outdir1 + '/filtered_contig_annotated.fasta'

gini_out = outdir1 + '/gini_index.txt'

file_vertices = outdir3 + '/Vertices_' + id + '.txt'
file_edges = outdir3 + '/Edges_' + id + '.txt'

file_vertices_comb = outdir3 + '/Vertices_comb_' + id + '.txt'
file_edges_comb = outdir3 + '/Edges_comb_' + id + '.txt'

layout_file = outdir3 + '/layout_'+ id + '.txt'

if(option == "1"):
    print('\nProcessing', '\t' + id + '\n')
    annot = Add_barcode_prefix(cellranger_annotation, annotation_, id)
    Format_fasta(cellranger_fasta, annot, annotated_fasta, id)
    RunIgBlastn(annotated_fasta, igdb)
    ParseIgBlast(igblast_out, annotated_fasta, germline, annotation_)
    changeo = Transfer_isotype(parsed_tab)
    changeo, gex_dict = Check_barcodes(changeo, sampleInfo_file, regex_pattern = '_[ATCG][ATCG]', barcode_columns = ('SANGER SAMPLE ID', 'GEX_SAMPLE_ID'))
    changeo = Check_contigs(changeo, obs = obs_file, obs_columns = ('index', 'filter'))
    changeo_filtered = Filter_changeo(changeo, filtered_tab, ambiguous_tab)

elapsed_time_secs = time() - dandelion_start_time
msg = "Execution took: %s secs (Wall clock time)\n" % timedelta(seconds=round(elapsed_time_secs))
print(msg)