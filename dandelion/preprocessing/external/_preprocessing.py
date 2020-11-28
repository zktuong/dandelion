#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 17:56:02
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-11-28 20:54:42

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
    reannotate with IgBLASTn

    Parameters
    ----------
    fasta
        fasta file
    igblast_db
        path to igblast database
    org
        organism
    loci
        ig or tr
    *args
        any arguments for ``AssignGenes.py``

    Returns
    -------
        igblast annotated file

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
    parses IgBLAST output to change-o format

    Parameters
    ----------
    fasta
        fasta file
    igblast_output
        igblast output file
    germline
        path to germline database
    fileformat
        format of output file
    org
        organism.
    verbose
        whether or not to print the files
    Returns
    -------
        change-o object
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

def tigger_genotype(data, germline=None, outdir=None, org = 'human', fileformat = 'airr', novel_ = 'NO', verbose = False):
    """
    reassignAlleles with TIgGER in R.

    Parameters
    ----------
    data
        Tabulated data, in Change-O (TAB) or AIRR (TSV) format.
    germline : PathLike, optional
        FASTA file containing IMGT-gapped V segment reference germlines.
        Defaults to $GERMLINE.
    outdir : PathLike,  optional
        Output directory. Will be created if it does not exist.
        Defaults to the current working directory.
    org : str
        organsim.
    fileformat
        Format for running tigger. Default is 'airr'. Also accepts 'changeo'.
    novel : str
        Whether or not to run novel allele discovery. Default is 'NO'.
    verbose : bool
        Whether or not to print the command used in the terminal. Default is False.
    Returns
    -------

    """

    start_time = time()
    env = os.environ.copy()
    if germline is None:
        try:
            gml = env['GERMLINE']
        except:
            raise OSError('Environmental variable GERMLINE must be set. Otherwise, please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files.')
        gml = gml+'imgt/'+org+'/vdj/imgt_'+org+'_IGHV.fasta'
    else:
        if os.path.isdir(germline):
            gml = germline.rstrip('/') + 'imgt_'+org+'_IGHV.fasta'
            if not os.path.isfile(gml):
                raise OSError("Input for germline is incorrect. Please rename IGHV germline file to '{}'. Otherwise, please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.".format(gml))
        elif type(germline) is not list:
            germline_ = [germline]
            if len(germline_) < 3:
                raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
            else:
                for x in germline_:
                    if not x.endswith('.fasta'):
                        raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                    if (os.path.isfile(x)) & ('ighv' in x.lower()):
                        gml = x
        else:
            if len(germline) < 3:
                raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
            else:
                for x in germline:
                    if not x.endswith('.fasta'):
                        raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                    if (os.path.isfile(x)) & ('ighv' in x.lower()):
                        gml = x

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

def recipe_scanpy_qc(self, max_genes=2500, min_genes=200, mito_cutoff=5, pval_cutoff=0.1, min_counts=None, max_counts=None, batch_term=None, blacklist=None):
    """
    Recipe for running a standard scanpy QC worflow.

    Parameters
    ----------
    adata : AnnData
        The (annotated) data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    max_genes : int
        Maximum number of genes expressed required for a cell to pass filtering. Default is 2500.
    min_genes : int
        Minimum number of genes expressed  required for a cell to pass filtering. Default is 200.
    mito_cutoff : float
        Maximum percentage mitochondrial content allowed for a cell to pass filtering. Default is 5.
    pval_cutoff : float
        Maximum Benjamini-Hochberg corrected p value from doublet detection protocol allowed for a cell to pass filtering. Default is 0.05.
    min_counts : int, optional
        Minimum number of counts required for a cell to pass filtering. Default is None.
    max_counts : int, optional
        Maximum number of counts required for a cell to pass filtering. Default is None.
    batch_term : str, optional
        If provided, will use sc.external.pp.bbknn for neighborhood construction.
    blacklist : sequence, optional
        If provided, will exclude these genes from highly variable genes list.
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
        if re.search('^TR[AB][VDJ]|^IG[HKL][VDJ]', i):
            _adata.var.at[i, 'highly_variable'] = False
        if blacklist is not None:
            if i in blacklist:
                _adata.var.at[i, 'highly_variable'] = False
    _adata = _adata[:, _adata.var['highly_variable']]
    sc.pp.scale(_adata, max_value=10)
    sc.tl.pca(_adata, svd_solver='arpack')
    if batch_term is None:
        sc.pp.neighbors(_adata, n_neighbors=10, n_pcs=50)
    else:
        sc.external.pp.bbknn(_adata, batch_key=batch_term)
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
    # let's do a one-sided test. the Bertie write-up does not address this but it makes sense
    pvals = 1-scipy.stats.norm.cdf(_adata.obs['scrublet_cluster_score'], loc=med, scale=1.4826*mad)
    _adata.obs['scrublet_score_bh_pval'] = bh(pvals)
    # threshold the p-values to get doublet calls.
    _adata.obs['is_doublet'] = _adata.obs['scrublet_score_bh_pval'] < pval_cutoff
    _adata.obs['is_doublet'] = _adata.obs['is_doublet'].astype('category')
    _adata.obs['filter_rna'] = (pd.Series([min_genes < n > max_genes for n in _adata.obs['n_genes']], index = _adata.obs.index)) | \
        (_adata.obs['percent_mito'] >= mito_cutoff) | \
            (_adata.obs['is_doublet'] == True)

    # removing columns that probably don't need anymore
    _adata.obs = _adata.obs.drop(['leiden', 'leiden_R', 'scrublet_cluster_score'], axis = 1)
    self.obs = _adata.obs.copy()