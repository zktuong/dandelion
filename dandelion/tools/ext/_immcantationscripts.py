#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-06-01 22:31:16
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-06-01 22:37:59

import os
import sys
from subprocess import run
import pandas as pd
import multiprocessing
from ...utilities._misc import *
from subprocess import run
from changeo.Gene import getGene

from rpy2.robjects.packages import importr, data
from rpy2.rinterface import NULL
from rpy2.robjects import pandas2ri
import rpy2.robjects
from plotnine import ggplot, geom_point, options, annotate, aes, xlab, ylab, facet_grid, theme_bw, geom_histogram, geom_vline, facet_grid, theme
import warnings
import numpy as np

def calculate_threshold(self, model=None, normalize_method=None, threshold_method=None, edge=None, cross=None, subsample=None, threshold_model=None, cutoff=None, sensitivity=None, specificity=None, ncpu=None, plot=True, plot_group=None, manual_threshold=None, figsize=(4.5, 2.5), *args):
    '''
    Wrappers for shazam's functions for calculating nearest neighbor distances for tuning clonal assignment.
        
    distToNearest
        Get non-zero distance of every heavy chain (IGH) sequence (as defined by sequenceColumn) to its nearest sequence in a partition of heavy chains sharing the same V gene, J gene, and junction length (VJL), or in a partition of single cells with heavy chains sharing the same heavy chain VJL combination, or of single cells with heavy and light chains sharing the same heavy chain VJL and light chain VJL combinations.
    findThreshold
        automtically determines an optimal threshold for clonal assignment of Ig sequences using a vector of nearest neighbor distances. It provides two alternative methods using either a Gamma/Gaussian Mixture Model fit (threshold_method="gmm") or kernel density fit (threshold_method="density").
    
    Parameters
    ----------
    model 
        underlying SHM model, which must be one of c("ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "hs1f_compat", "m1n_compat").
    normalize_method  
        method of normalization. The default is "len", which divides the distance by the length of the sequence group. If "none" then no normalization if performed.
    threshold_method
        string defining the method to use for determining the optimal threshold. One of "gmm" or "density".
    edge
        upper range as a fraction of the data density to rule initialization of Gaussian fit parameters. Default value is 0.9 (or 90). Applies only when threshold_method="density".
    cross
        supplementary nearest neighbor distance vector output from distToNearest for initialization of the Gaussian fit parameters. Applies only when method="gmm".
    subsample
        maximum number of distances to subsample to before threshold detection.
    threshold_model
        allows the user to choose among four possible combinations of fitting curves: "norm-norm", "norm-gamma", "gamma-norm", and "gamma-gamma". Applies only when method="gmm".
    cutoff
        method to use for threshold selection: the optimal threshold "opt", the intersection point of the two fitted curves "intersect", or a value defined by user for one of the sensitivity or specificity "user". Applies only when method="gmm".
    sensitivity
        sensitivity required. Applies only when method="gmm" and cutoff="user".
    specificity
        specificity required. Applies only when method="gmm" and cutoff="user".
    *args
        passed to shazam's distToNearest function
    '''

    sh = importr('shazam')
    dat = load_data(self)
    warnings.filterwarnings("ignore")

    if 'v_call_genotyped' in dat.columns:
        v_call = 'v_call_genotyped'
    else:
        v_call = 'v_call'
    
    if model is None:
        model_ = 'ham'
    else:
        model_ = model
    
    if normalize_method is None:
        norm_ = 'len'
    else:
        norm_ = normalize_method
    
    if threshold_method is None:
        threshold_method_ = "density"
    else:
        threshold_method_ = threshold_method
    
    if subsample is None:
        subsample_ = NULL
    else:
        subsample_ = subsample
    
    if ncpu is None:
        ncpu_ = multiprocessing.cpu_count()
    else:
        ncpu_ = ncpu

    dat_h = dat[dat['locus'] == 'IGH']
    
    try:
        dat_h_r = pandas2ri.py2rpy(dat_h)
    except:
        dat_h = dat_h.fillna('')
        dat_h_r = pandas2ri.py2rpy(dat_h)

    dist_ham = sh.distToNearest(dat_h_r, vCallColumn=v_call, model=model_, normalize=norm_, nproc=ncpu_, *args)
    # Find threshold using density method
    c = rpy2.robjects.StrVector(['dist_nearest'])

    if threshold_method_ is 'density':
        if edge is None:
            edge_ = 0.9
        else:
            edge_ = edge
            
        dist_threshold = sh.findThreshold(dist_ham.rx(True, c), method=threshold_method_, subsample = subsample_, edge = edge_)
    else:
        if threshold_model is None:
            threshold_model_ = "gamma-gamma"
        else:
            threshold_model_ = threshold_model
        
        if cross is None:
            cross_ = NULL
        else:
            cross_ = cross
        
        if cutoff is None:
            cutoff_ = 'optimal'
        else:
            cutoff_ = cutoff
        
        if sensitivity is None:
            sen_ = NULL
        else:
            sen_ = sensitivity
        
        if specificity is None:
            spc_ = NULL
        else:
            spc_ = specificity        
        dist_threshold = sh.findThreshold(dist_ham.rx(True, c), method=threshold_method, model = threshold_model_, cross = cross_, subsample = subsample_, cutoff = cutoff_, sen = sen_, spc = spc_)        

    threshold=np.array(dist_threshold.slots['threshold'])[0]
    
    dist_ham = pandas2ri.rpy2py_dataframe(dist_ham)
    
    if plot:
        options.figure_size = figsize
        if plot_group is None:
            plot_group = 'sample_id'
        else:
            plot_group = plot_group
        if manual_threshold is None:
            tr = threshold
        else:
            tr = manual_threshold            
        p = (ggplot(dist_ham, aes('dist_nearest', fill=str(plot_group)))
             + theme_bw() 
             + xlab("Grouped Hamming distance")
             + ylab("Count")
             + geom_histogram(binwidth = 0.01)
             + geom_vline(xintercept = tr, linetype = "dashed", color="blue", size=0.5)
             + annotate('text', x=tr+0.02, y = 10, label='Threshold:\n' + str(np.around(tr, decimals=2)), size = 8, color = 'Blue', hjust = 'left')
             + facet_grid('~'+str(plot_group), scales="free_y")
             + theme(legend_position = 'none'))        
        return(p)            
    else:
        print('Automatic Threshold : '+str(np.around(threshold, decimals=2), '\n method = '+str(threshold_method)))

def define_clones(data, dist, action = 'set', model = 'ham', norm = 'len', doublets='drop', fileformat='airr', ncpu = None, dirs = None, verbose = False, *args):
    """    
    """
    if ncpu is None:
        nproc=multiprocessing.cpu_count()
    else:
        nproc=ncpu

    dat = load_data(data)
    dat_h = dat[dat['locus'] == 'IGH']
    dat_l = dat[dat['locus'].isin(['IGK', 'IGL'])]
    
    if dirs is None:
        tmpFolder = "{}/tmp".format(os.path.dirname(data))
        outFolder = "{}".format(os.path.dirname(data))
    else:
        tmpFolder = dirs
        outFolder = dirs

    if not os.path.exists(tmpFolder):
        os.makedirs(tmpFolder)

    h_file1 = "{}/{}_heavy-clone.tsv".format(tmpFolder, os.path.basename(data).split('.tsv')[0])
    h_file2 = "{}/{}_heavy-clone.tsv".format(tmpFolder, os.path.basename(data).split('.tsv')[0])
    l_file = "{}/{}_light.tsv".format(tmpFolder, os.path.basename(data).split('.tsv')[0])
    outfile = "{}/{}_clone.tsv".format(outFolder, os.path.basename(data).split('.tsv')[0])

    dat_h.to_csv(h_file1, sep = '\t', index = False)
    dat_l.to_csv(l_file, sep = '\t', index = False)
    
    if 'v_call_genotyped' in dat.columns:
        v_field = 'v_call_genotyped'
    else:
        v_field = 'v_call'

    cmd = ['DefineClones.py',
            '-d', h_file1,
            '-o', h_file2,
            '--act', action,
            '--model', model,
            '--norm', norm,
            '--dist', str(dist),
            '--nproc', str(nproc),
            '--vf', v_field,
            *args]    
    
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

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd)
    lightCluster(h_file2, l_file, outfile, doublets=doublets, fileformat=fileformat)
