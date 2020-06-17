#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 17:56:02
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-06-17 13:42:50

import os
from subprocess import run
from datetime import timedelta
from time import time
from collections import OrderedDict
from time import time

# Presto and changeo imports
from presto.IO import printLog, printProgress, printError, printWarning
from changeo.Alignment import gapV
from changeo.Defaults import default_format, default_out_args
from changeo.Gene import buildGermline
from changeo.IO import countDbFile, getFormatOperators, getOutputHandle, readGermlines, checkFields

def assigngenes_igblast(fasta, igblast_db = None, org = 'human', loci = 'ig', fileformat = 'airr', verbose = False, outputfolder = None):
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
    fileformat: str (Default: 'airr')
        format of the output data. Default is 'airr'. changeo' will trigger legacy changeo format.
    *args
        any arguments for ``AssignGenes.py``

    Returns
    -------
        igblast annotated file

    """

    file_format_dict = {'changeo':'blast', 'airr':'airr'}
    env = os.environ.copy()
    if igblast_db is None:
        try:
            igdb = env['IGDATA']
        except:
            raise OSError('Environmental variable IGDATA must be set. Otherwise, please provide path to igblast database')
    else:
        env['IGDATA'] = igblast_db
        igdb = env['IGDATA']

    cmd = ['AssignGenes.py', 'igblast',
           '-s', fasta,
           '-b', igdb,
           '--organism', org,
           '--loci', loci,
           '--format', file_format_dict[fileformat]
           ]

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd, env=env) # logs are printed to terminal

    informat_dict = {'changeo':'_igblast.fmt7'}

    if fileformat == 'changeo':
        if outputfolder is None:
            outfolder = os.path.dirname(fasta)+'/tmp'
        else:
            outfolder = str(outputfolder)
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        in_file = fasta.split('.fasta')[0] + informat_dict[fileformat]
        outfile = os.path.basename(fasta).split('.fasta')[0] + informat_dict[fileformat]
        out_file = "{}/{}".format(outfolder, outfile)
        os.replace(in_file, out_file)

def makedb_igblast(fasta, igblast_output = None, germline = None, org = 'human', outputfolder = None, extended = False, verbose = False):
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
    outputfolder
        if specified, the location of output files.
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
            raise OSError('Environmental variable GERMLINE must be set. Otherwise, please provide path to germline fasta files')
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
            '--asis-id',
            '--extended']
    else:
        cmd = ['MakeDb.py', 'igblast',
               '-i', igbo,
               '-s', fasta,
               '-r', gml,
               '--10x', cellranger_annotation,
               '--asis-id']

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd, env=env) # logs are printed to terminal

    new_file = igbo.replace('.fmt7', '_db-pass.tsv')
    outfolder = os.path.dirname(new_file).split('/tmp')[0]
    out_file = "{}/{}".format(outfolder, os.path.basename(new_file))
    os.replace(new_file, out_file)

def tigger_genotype(data, germline=None, outdir=None, org = 'human', fileformat = 'airr', verbose = False):
    """
    reassignAlleles with TIgGER in R.

    Parameters
    ----------
    data
        Tabulated data, in Change-O (TAB) or AIRR (TSV) format.
    germline
        FASTA file containing IMGT-gapped V segment reference germlines.
        Defaults to $GERMLINE.
    outdir
        Output directory. Will be created if it does not exist.
        Defaults to the current working directory.
    *args
        any arguments for ``tigger-genotype.R``.

    Returns
    -------

    """

    start_time = time()
    env = os.environ.copy()
    if germline is None:
        try:
            gml = env['GERMLINE']
        except:
            raise OSError('Environmental variable GERMLINE must be set. Otherwise, please provide path to germline fasta files.')
        gml = gml+'imgt/'+org+'/vdj/imgt_'+org+'_IGHV.fasta'
    else:
        env['GERMLINE'] = germline
        gml = germline

    if not gml.endswith('.fasta'):
        raise OSError('Input for germline is incorrect. Please provide path to germline IGHV fasta file.')

    if outdir is not None:
        out_dir = outdir + '/'
    else:
        out_dir = os.path.dirname(data)

    cmd = ['tigger-genotype.R',
           '-d', data,
           '-r', gml,
           '-n', os.path.basename(data).split('.tsv')[0],
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

def insertGaps(db_file, references=None, format=default_format,
               out_file=None, out_args=default_out_args):
    """
    Inserts IMGT numbering into V fields

    Arguments:
      db_file : the database file name.
      references : folder with germline repertoire files. If None, do not updated alignment columns wtih IMGT gaps.
      format : input format.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
     str : output file name
    """
    log = OrderedDict()
    log['START'] = 'insertGaps'
    log['COMMAND'] = 'insertGaps'
    log['FILE'] = os.path.basename(db_file)
    # printLog(log)

    # Define format operators
    try:
        reader, writer, schema = getFormatOperators(format)
    except ValueError:
        printError('Invalid format %s.' % format)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = reader(db_handle)

    # Check for required columns
    try:
        required = ['sequence_imgt', 'v_germ_start_imgt']
        checkFields(required, db_iter.fields, schema=schema)
    except LookupError as e:
        printError(e)

    # Load references
    reference_dict = readGermlines(references)

    # Check for IMGT-gaps in germlines
    if all('...' not in x for x in reference_dict.values()):
        printWarning('Germline reference sequences do not appear to contain IMGT-numbering spacers. Results may be incorrect.')

    # Open output writer
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='gap', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=schema.out_type)
    pass_writer = writer(pass_handle, fields=db_iter.fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count = pass_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        # printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1
        # Update IMGT fields
        imgt_dict = correctIMGTFields(rec, reference_dict)
        # Write records
        if imgt_dict is not None:
            pass_count += 1
            rec.setDict(imgt_dict, parse=False)
            pass_writer.writeReceptor(rec)

    # Print counts
    # printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['PASS'] = pass_count
    log['FAIL'] = rec_count - pass_count
    log['END'] = 'insertGaps'
    # printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name

def correctIMGTFields(receptor, references):
    """
    Add IMGT-gaps to IMGT fields in a Receptor object

    Arguments:
      receptor (changeo.Receptor.Receptor): Receptor object to modify.
      references (dict): dictionary of IMGT-gapped references sequences.

    Returns:
      changeo.Receptor.Receptor: modified Receptor with IMGT-gapped fields.
    """
    # Initialize update object
    imgt_dict = {'sequence_imgt': None,
                 'v_germ_start_imgt': None,
                 'v_germ_length_imgt': None,
                 'germline_imgt': None}

    try:
        if not all([receptor.sequence_imgt,
                    receptor.v_germ_start_imgt,
                    receptor.v_germ_length_imgt,
                    receptor.v_call]):
            raise AttributeError
    except AttributeError:
        return None

    # Update IMGT fields
    try:
        gapped = gapV(receptor.sequence_imgt,
                      receptor.v_germ_start_imgt,
                      receptor.v_germ_length_imgt,
                      receptor.v_call,
                      references)
    except KeyError as e:
        printWarning(e)
        return None

    # Verify IMGT-gapped sequence and junction concur
    try:
        check = (receptor.junction == gapped['sequence_imgt'][309:(309 + receptor.junction_length)])
    except TypeError:
        check = False
    if not check:
        return None

    # Rebuild germline sequence
    __, germlines, __ = buildGermline(receptor, references)
    if germlines is None:
        return None
    else:
        gapped['germline_imgt'] = germlines['full']

    # Update return object
    imgt_dict.update(gapped)

    return imgt_dict