#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 17:56:02
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-12-29 20:58:23

import os
import pandas as pd
import numpy as np
from subprocess import run
from datetime import timedelta
from time import time
from collections import OrderedDict
from time import time
from ...utilities._utilities import *
import scanpy as sc
import scipy.stats
import re

def assigngenes_igblast(fasta, igblast_db = None, org = 'human', loci = 'ig', verbose = False):
    """
    Reannotate with IgBLASTn.

    Parameters
    ----------
    fasta : PathLike
        fasta file for reannotation.
    igblast_db : PathLike, optional
        path to igblast database.
    org : str
        organism for germline sequences.
    loci : str
        `ig` or `tr` mode for running igblastn.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.
    Returns
    -------
    """
    env = os.environ.copy()
    if igblast_db is None:
        try:
            igdb = env['IGDATA']
        except:
            raise OSError('Environmental variable IGDATA must be set. Otherwise, please provide path to igblast database')
    else:
        env['IGDATA'] = igblast_db
        igdb = env['IGDATA']

    outfolder = os.path.abspath(os.path.dirname(fasta))+'/tmp'
    informat_dict = {'blast':'_igblast.fmt7', 'airr':'_igblast.tsv'}
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    
    for fileformat in ['blast', 'airr']:
        outfile = os.path.basename(fasta).split('.fasta')[0] + informat_dict[fileformat]
        cmd = ['AssignGenes.py', 'igblast',
               '-s', fasta,
               '-b', igdb,
               '--organism', org,
               '--loci', loci,
               '--format', fileformat,
               '-o', "{}/{}".format(outfolder, outfile)
           ]

        if verbose:
            print('Running command: %s\n' % (' '.join(cmd)))
        run(cmd, env=env) # logs are printed to terminal

def makedb_igblast(fasta, igblast_output = None, germline = None, org = 'human', extended = True, verbose = False):
    """
    Parses IgBLAST output to airr format.

    Parameters
    ----------
    fasta : PathLike
        fasta file use for reannotation.
    igblast_output : PathLike, optional
        igblast output file.
    germline : PathLike, optional
        path to germline database.    
    org : str
        organism of germline sequences.
    extended : bool
        whether or not to parse extended 10x annotations. Default is True.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.
    Returns
    -------
    """
    env = os.environ.copy()
    if germline is None:
        try:
            gml = env['GERMLINE']
        except:
            raise OSError('Environmental variable GERMLINE must be set. Otherwise, please provide path to folder containing germline fasta files.')
        gml = gml+'imgt/'+org+'/vdj/'
    else:
        env['GERMLINE'] = germline
        gml = germline

    if igblast_output is None:
        indir = os.path.dirname(fasta)+'/tmp'
        infile = os.path.basename(fasta).split('.fasta')[0] + '_igblast.fmt7'
        igbo = "{}/{}".format(indir, infile)
    else:
        igbo = igblast_output

    cellranger_annotation = "{}/{}".format(os.path.dirname(fasta), os.path.basename(fasta).replace('.fasta', '_annotations.csv'))

    if extended:
        cmd = ['MakeDb.py', 'igblast',
            '-i', igbo,
            '-s', fasta,
            '-r', gml,
            '--10x', cellranger_annotation,
            '--extended']
    else:
        cmd = ['MakeDb.py', 'igblast',
               '-i', igbo,
               '-s', fasta,
               '-r', gml,
               '--10x', cellranger_annotation]

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd, env=env) # logs are printed to terminal

def parsedb_heavy(db_file, verbose = False):
    """
    Parses AIRR table (heavy chain contigs only).

    Parameters
    ----------
    db_file : PathLike
        path to AIRR table.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.
    Returns
    -------    
    """    
    outname = os.path.basename(db_file).split('.tsv')[0] + '_heavy'
    
    cmd = ['ParseDb.py', 'select',
           '-d', db_file,
           '-f', 'locus',
           '-u', 'IGH',
           '--logic', 'all',
           '--regex',
           '--outname', outname 
           ]
     
    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd) # logs are printed to terminal

def parsedb_light(db_file, verbose = False):
    """
    Parses AIRR table (light chain contigs only).

    Parameters
    ----------
    db_file : PathLike
        path to AIRR table.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.
    Returns
    -------    
    """    
    outname = os.path.basename(db_file).split('.tsv')[0] + '_light'
    
    cmd = ['ParseDb.py', 'select',
           '-d', db_file,
           '-f', 'locus',
           '-u', 'IG[LK]',
           '--logic', 'all',
           '--regex',
           '--outname', outname 
           ]
     
    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd) # logs are printed to terminal

def creategermlines(db_file, germtypes = None, germline = None, org = 'human', genotype_fasta = None, v_field = None, cloned = False, mode = None, verbose = False):
    """
    Wrapper for CreateGermlines.py for reconstructing germline sequences,

    Parameters
    ----------
    db_file : PathLike
        path to AIRR table.
    germtypes : str, optional
        germline type for reconstuction.
    germline : PathLike, optional
        location to germline fasta files.
    org : str
        organism for germline sequences.
    genotype_fasta : PathLike, optional
        location to corrected v germine fasta file.
    v_field : str, optional
        name of column for v segment to perform reconstruction.
    cloned : bool
        whether or not to run with cloned option.
    mode : str, optional
        whether to run on heavy or light mode. If left as None, heavy and light will be run together.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.
    Returns
    -------    
    """
    env = os.environ.copy()
    if germline is None:
        try:
            gml = env['GERMLINE']
        except:
            raise OSError('Environmental variable GERMLINE must be set. Otherwise, please provide path to folder containing germline fasta files.')
        gml = gml+'imgt/'+org+'/vdj/'
    else:
        env['GERMLINE'] = germline
        gml = germline
    
    if germtypes is None:
        germ_type = 'dmask'
    else:
        germ_type = germtypes

    if cloned:        
        if mode == 'heavy':
            print('            Reconstructing heavy chain {} germline sequences with {} for each clone.'.format(germ_type, v_field))
            if genotype_fasta is None:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file, 
                        '-g', germ_type, 
                        '--cloned',
                        '-r', gml+'/imgt_'+org+'_IGHV.fasta', gml+'/imgt_'+org+'_IGHD.fasta', gml+'/imgt_'+org+'_IGHJ.fasta',
                        '--vf', v_field
                        ]
                else:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '--cloned',
                        '-r', gml,
                        '--vf', v_field
                        ]
            else:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '--cloned',
                        '-r', genotype_fasta, gml+'/imgt_'+org+'_IGHD.fasta', gml+'/imgt_'+org+'_IGHJ.fasta',
                        '--vf', v_field
                        ]
                else:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '--cloned',
                        '-r', genotype_fasta, gml,
                        '--vf', v_field
                        ]
        elif mode == 'light':
            print('            Reconstructing light chain {} germline sequences with {} for each clone.'.format(germ_type, v_field))
            if germline is None:
                cmd = ['CreateGermlines.py',
                    '-d', db_file,
                    '-g', germ_type,
                    '--cloned',
                    '-r', gml+'/imgt_'+org+'_IGKV.fasta', gml+'/imgt_'+org+'_IGKJ.fasta', gml+'/imgt_'+org+'_IGLV.fasta', gml+'/imgt_'+org+'_IGLJ.fasta',
                    '--vf', v_field
                    ]
            else:
                cmd = ['CreateGermlines.py',
                    '-d', db_file,
                    '-g', germ_type,
                    '--cloned',
                    '-r', gml,
                    '--vf', v_field
                    ]
        elif mode is None:
            print('            Reconstructing {} germline sequences with {} for each clone.'.format(germ_type, v_field))
            if genotype_fasta is None:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '--cloned',
                        '-r', gml+'/imgt_'+org+'_IGHV.fasta', gml+'/imgt_'+org+'_IGHD.fasta', gml+'/imgt_'+org+'_IGHJ.fasta', gml+'/imgt_'+org+'_IGKV.fasta', gml+'/imgt_'+org+'_IGKJ.fasta', gml+'/imgt_'+org+'_IGLV.fasta', gml+'/imgt_'+org+'_IGLJ.fasta',
                        '--vf', v_field
                        ]
                else:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '--cloned',
                        '-r', gml,
                        '--vf', v_field
                        ]
            else:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '--cloned',
                        '-r', genotype_fasta, gml+'/imgt_'+org+'_IGHD.fasta', gml+'/imgt_'+org+'_IGHJ.fasta', gml+'/imgt_'+org+'_IGKV.fasta', gml+'/imgt_'+org+'_IGKJ.fasta', gml+'/imgt_'+org+'_IGLV.fasta', gml+'/imgt_'+org+'_IGLJ.fasta',
                        '--vf', v_field
                        ]
                else:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '--cloned',
                        '-r', genotype_fasta, gml,
                        '--vf', v_field
                        ]
    else:        
        if mode == 'heavy':
            print('            Reconstructing heavy chain {} germline sequences with {}.'.format(germ_type, v_field))
            if genotype_fasta is None:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '-r', gml+'/imgt_'+org+'_IGHV.fasta', gml+'/imgt_'+org+'_IGHD.fasta', gml+'/imgt_'+org+'_IGHJ.fasta',
                        '--vf', v_field
                        ]
                else:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '-r', gml,
                        '--vf', v_field
                        ]
            else:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '-r', genotype_fasta, gml+'/imgt_'+org+'_IGHD.fasta', gml+'/imgt_'+org+'_IGHJ.fasta',
                        '--vf', v_field
                        ]
                else:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '-r', genotype_fasta, gml,
                        '--vf', v_field
                        ]
        elif mode == 'light':
            print('            Reconstructing light chain {} germline sequences with {}.'.format(germ_type, v_field))
            if germline is None:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '-r', gml+'/imgt_'+org+'_IGKV.fasta', gml+'/imgt_'+org+'_IGKJ.fasta', gml+'/imgt_'+org+'_IGLV.fasta', gml+'/imgt_'+org+'_IGLJ.fasta',
                        '--vf', v_field
                        ]
            else:
                cmd = ['CreateGermlines.py',
                    '-d', db_file,
                    '-g', germ_type,
                    '-r', gml,
                    '--vf', v_field
                    ]
        elif mode is None:
            print('            Reconstructing {} germline sequences with {} for each clone.'.format(germ_type, v_field))
            if genotype_fasta is None:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '-r', gml+'/imgt_'+org+'_IGHV.fasta', gml+'/imgt_'+org+'_IGHD.fasta', gml+'/imgt_'+org+'_IGHJ.fasta', gml+'/imgt_'+org+'_IGKV.fasta', gml+'/imgt_'+org+'_IGKJ.fasta', gml+'/imgt_'+org+'_IGLV.fasta', gml+'/imgt_'+org+'_IGLJ.fasta',
                        '--vf', v_field
                        ]
                else:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '-r', gml,
                        '--vf', v_field
                        ]
            else:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '-r', genotype_fasta, gml+'/imgt_'+org+'_IGHD.fasta', gml+'/imgt_'+org+'_IGHJ.fasta', gml+'/imgt_'+org+'_IGKV.fasta', gml+'/imgt_'+org+'_IGKJ.fasta', gml+'/imgt_'+org+'_IGLV.fasta', gml+'/imgt_'+org+'_IGLJ.fasta',
                        '--vf', v_field
                        ]
                else:
                    cmd = ['CreateGermlines.py',
                        '-d', db_file,
                        '-g', germ_type,
                        '-r', genotype_fasta, gml,
                        '--vf', v_field
                        ]

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd, env = env) # logs are printed to terminal

def tigger_genotype(data, v_germline=None, outdir=None, org = 'human', fileformat = 'airr', novel_ = 'YES', verbose = False):
    """
    Reassign alleles with TIgGER in R.

    Parameters
    ----------
    data : PathLike
        vdj tabulated data, in Change-O (TAB) or AIRR (TSV) format.
    germline : PathLike, optional
        fasta file containing IMGT-gapped V segment reference germlines. Defaults to $GERMLINE.
    outdir : PathLike,  optional
        output directory. Will be created if it does not exist. Defaults to the current working directory.
    org : str
        organism for germline sequences.
    fileformat : str
        format for running tigger. Default is 'airr'. Also accepts 'changeo'.
    novel : str
        whether or not to run novel allele discovery. Default is 'YES'.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.
    Returns
    -------
    """
    start_time = time()
    env = os.environ.copy()
    if v_germline is None:
        try:
            gml = env['GERMLINE']
        except:
            raise OSError('Environmental variable GERMLINE is not set. Please provide either the path to the folder containing the germline IGHV fasta file, or direct path to the germline IGHV fasta file.')
        gml = gml+'imgt/'+org+'/vdj/imgt_'+org+'_IGHV.fasta'
    else:
        if os.path.isdir(v_germline):
            gml = v_germline.rstrip('/') + 'imgt_'+org+'_IGHV.fasta'
            if not os.path.isfile(gml):
                raise OSError("Input for germline is incorrect. Please rename IGHV germline file to '{}'. Otherwise, please provide path to folder containing the germline IGHV fasta file, or direct path to the germline IGHV fasta file.".format(gml))        
        else:            
            if not v_germline.endswith('.fasta'):
                raise OSError('Input for germline is incorrect {}. Please provide path to folder containing the germline IGHV fasta file, or direct path to the germline IGHV fasta file.'.format(v_germline))
            if (os.path.isfile(v_germline)) & ('ighv' in v_germline.lower()):
                gml = v_germline

    if outdir is not None:
        out_dir = outdir + '/'
    else:
        out_dir = os.path.dirname(data)

    cmd = ['tigger-genotype.R',
           '-d', data,
           '-r', gml,
           '-n', os.path.basename(data).split('.tsv')[0],
           '-N', novel_,
           '-o', out_dir,
           '-f', fileformat]

    print('      Reassigning alleles')
    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd, env=env) # logs are printed to terminal
    elapsed_time_secs = time() - start_time
    msg = "tigger-genotype execution took: %s secs (Wall clock time)\n" % timedelta(seconds=round(elapsed_time_secs))
    if verbose:
        print(msg)

## commented out originally in ConvertDb, not sure if it works properly
# def insertGaps(db_file, references=None, format=default_format,
#                out_file=None, out_args=default_out_args):
#     """
#     Inserts IMGT numbering into V fields

#     Arguments:
#       db_file : the database file name.
#       references : folder with germline repertoire files. If None, do not updated alignment columns wtih IMGT gaps.
#       format : input format.
#       out_file : output file name. Automatically generated from the input file if None.
#       out_args : common output argument dictionary from parseCommonArgs.

#     Returns:
#      str : output file name
#     """
#     log = OrderedDict()
#     log['START'] = 'insertGaps'
#     log['COMMAND'] = 'insertGaps'
#     log['FILE'] = os.path.basename(db_file)
#     # printLog(log)

#     # Define format operators
#     try:
#         reader, writer, schema = getFormatOperators(format)
#     except ValueError:
#         printError('Invalid format %s.' % format)

#     # Open input
#     db_handle = open(db_file, 'rt')
#     db_iter = reader(db_handle)

#     # Check for required columns
#     try:
#         required = ['sequence_imgt', 'v_germ_start_imgt']
#         checkFields(required, db_iter.fields, schema=schema)
#     except LookupError as e:
#         printError(e)

#     # Load references
#     reference_dict = readGermlines(references)

#     # Check for IMGT-gaps in germlines
#     if all('...' not in x for x in reference_dict.values()):
#         printWarning('Germline reference sequences do not appear to contain IMGT-numbering spacers. Results may be incorrect.')

#     # Open output writer
#     if out_file is not None:
#         pass_handle = open(out_file, 'w')
#     else:
#         pass_handle = getOutputHandle(db_file, out_label='gap', out_dir=out_args['out_dir'],
#                                       out_name=out_args['out_name'], out_type=schema.out_type)
#     pass_writer = writer(pass_handle, fields=db_iter.fields)

#     # Count records
#     result_count = countDbFile(db_file)

#     # Iterate over records
#     start_time = time()
#     rec_count = pass_count = 0
#     for rec in db_iter:
#         # Print progress for previous iteration
#         # printProgress(rec_count, result_count, 0.05, start_time=start_time)
#         rec_count += 1
#         # Update IMGT fields
#         imgt_dict = correctIMGTFields(rec, reference_dict)
#         # Write records
#         if imgt_dict is not None:
#             pass_count += 1
#             rec.setDict(imgt_dict, parse=False)
#             pass_writer.writeReceptor(rec)

#     # Print counts
#     # printProgress(rec_count, result_count, 0.05, start_time=start_time)
#     log = OrderedDict()
#     log['OUTPUT'] = os.path.basename(pass_handle.name)
#     log['RECORDS'] = rec_count
#     log['PASS'] = pass_count
#     log['FAIL'] = rec_count - pass_count
#     log['END'] = 'insertGaps'
#     # printLog(log)

#     # Close file handles
#     pass_handle.close()
#     db_handle.close()

#     return pass_handle.name

# def correctIMGTFields(receptor, references):
#     """
#     Add IMGT-gaps to IMGT fields in a Receptor object

#     Arguments:
#       receptor (changeo.Receptor.Receptor): Receptor object to modify.
#       references (dict): dictionary of IMGT-gapped references sequences.

#     Returns:
#       changeo.Receptor.Receptor: modified Receptor with IMGT-gapped fields.
#     """
#     # Initialize update object
#     imgt_dict = {'sequence_imgt': None,
#                  'v_germ_start_imgt': None,
#                  'v_germ_length_imgt': None,
#                  'germline_imgt': None}

#     try:
#         if not all([receptor.sequence_imgt,
#                     receptor.v_germ_start_imgt,
#                     receptor.v_germ_length_imgt,
#                     receptor.v_call]):
#             raise AttributeError
#     except AttributeError:
#         return None

#     # Update IMGT fields
#     try:
#         gapped = gapV(receptor.sequence_imgt,
#                       receptor.v_germ_start_imgt,
#                       receptor.v_germ_length_imgt,
#                       receptor.v_call,
#                       references)
#     except KeyError as e:
#         printWarning(e)
#         return None

#     # Verify IMGT-gapped sequence and junction concur
#     try:
#         check = (receptor.junction == gapped['sequence_imgt'][309:(309 + receptor.junction_length)])
#     except TypeError:
#         check = False
#     if not check:
#         return None

#     # Rebuild germline sequence
#     __, germlines, __ = buildGermline(receptor, references)
#     if germlines is None:
#         return None
#     else:
#         gapped['germline_imgt'] = germlines['full']

#     # Update return object
#     imgt_dict.update(gapped)

#     return imgt_dict

def recipe_scanpy_qc(self, max_genes=2500, min_genes=200, mito_cutoff=5, pval_cutoff=0.1, min_counts=None, max_counts=None, blacklist=None):
    """
    Recipe for running a standard scanpy QC workflow.

    Parameters
    ----------
    adata : AnnData
        annotated data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    max_genes : int
        naximum number of genes expressed required for a cell to pass filtering. Default is 2500.
    min_genes : int
        minimum number of genes expressed  required for a cell to pass filtering. Default is 200.
    mito_cutoff : float
        maximum percentage mitochondrial content allowed for a cell to pass filtering. Default is 5.
    pval_cutoff : float
        maximum Benjamini-Hochberg corrected p value from doublet detection protocol allowed for a cell to pass filtering. Default is 0.05.
    min_counts : int, optional
        minimum number of counts required for a cell to pass filtering. Default is None.
    max_counts : int, optional
        maximum number of counts required for a cell to pass filtering. Default is None.
    blacklist : sequence, optional
        if provided, will exclude these genes from highly variable genes list.
    Returns
    -------
        `AnnData` of shape n_obs × n_vars where obs now contain filtering information. Rows correspond to cells and columns to genes.
    """
    _adata = self.copy()
    # run scrublet    
    try:
        import scrublet as scr
    except:
        raise ImportError('Please install scrublet with pip install scrublet.')

    scrub = scr.Scrublet(_adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    _adata.obs['scrublet_score'] = doublet_scores
    # overcluster prep. run basic scanpy pipeline
    sc.pp.filter_cells(_adata, min_genes = 0)
    mito_genes = _adata.var_names.str.startswith('MT-')
    _adata.obs['percent_mito'] = np.sum(_adata[:, mito_genes].X, axis = 1) / np.sum(_adata.X, axis = 1)*100
    _adata.obs['n_counts'] = _adata.X.sum(axis = 1)
    sc.pp.normalize_total(_adata, target_sum = 1e4)
    sc.pp.log1p(_adata)
    sc.pp.highly_variable_genes(_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    for i in _adata.var.index:
        if re.search('^TR[AB][VDJ]|^IG[HKL][VDJC]', i):
            _adata.var.at[i, 'highly_variable'] = False
        if blacklist is not None:
            if i in blacklist:
                _adata.var.at[i, 'highly_variable'] = False
    _adata = _adata[:, _adata.var['highly_variable']]
    sc.pp.scale(_adata, max_value=10)
    sc.tl.pca(_adata, svd_solver='arpack')
    sc.pp.neighbors(_adata, n_neighbors=10, n_pcs=50)
    # overclustering proper - do basic clustering first, then cluster each cluster
    sc.tl.leiden(_adata)
    for clus in list(np.unique(_adata.obs['leiden']))[0]:
        sc.tl.leiden(_adata, restrict_to=('leiden',[clus]), key_added = 'leiden_R')
    for clus in list(np.unique(_adata.obs['leiden']))[1:]: # weird how the new anndata/scanpy is forcing this
        sc.tl.leiden(_adata, restrict_to=('leiden_R',[clus]), key_added = 'leiden_R')
    # compute the cluster scores - the median of Scrublet scores per overclustered cluster
    for clus in np.unique(_adata.obs['leiden_R']):
        _adata.obs.loc[_adata.obs['leiden_R']==clus, 'scrublet_cluster_score'] = \
            np.median(_adata.obs.loc[_adata.obs['leiden_R']==clus, 'scrublet_score'])
    # now compute doublet p-values. figure out the median and mad (from above-median values) for the distribution
    med = np.median(_adata.obs['scrublet_cluster_score'])
    mask = _adata.obs['scrublet_cluster_score']>med
    mad = np.median(_adata.obs['scrublet_cluster_score'][mask]-med)
    # 1 sided test for catching outliers
    pvals = 1-scipy.stats.norm.cdf(_adata.obs['scrublet_cluster_score'], loc=med, scale=1.4826*mad)
    _adata.obs['scrublet_score_bh_pval'] = bh(pvals)
    # threshold the p-values to get doublet calls.
    _adata.obs['is_doublet'] = _adata.obs['scrublet_score_bh_pval'] < pval_cutoff
    _adata.obs['is_doublet'] = _adata.obs['is_doublet'].astype('category')
    _adata.obs['filter_rna'] = (pd.Series([min_genes < n > max_genes for n in _adata.obs['n_genes']], index = _adata.obs.index)) | \
        (_adata.obs['percent_mito'] >= mito_cutoff) | \
            (_adata.obs['is_doublet'] == True)

    # removing columns that probably don't need anymore
    _adata.obs = _adata.obs.drop(['leiden', 'leiden_R', 'scrublet_cluster_score', 'scrublet_score_bh_pval'], axis = 1)
    self.obs = _adata.obs.copy()