#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 17:56:02
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-06-07 21:47:14

import sys
import os
import pandas as pd
from subprocess import run
from tqdm import tqdm
import multiprocessing
from joblib import Parallel, delayed
from collections import OrderedDict
from time import sleep
from ..utilities._misc import *
from .ext._immcantationscripts import assigngenes_igblast, makedb_igblast, tigger_genotype, insertGaps
from plotnine import ggplot, geom_bar, ggtitle, scale_fill_manual, coord_flip, options, element_blank, aes, xlab, ylab, facet_grid, theme_classic, theme
from changeo.Gene import buildGermline
from changeo.IO import countDbFile, getDbFields, getFormatOperators, readGermlines, checkFields
from changeo.Receptor import AIRRSchema, ChangeoSchema, Receptor, ReceptorData
import re
import scanpy as sc
import numpy as np
import scipy.stats
import scrublet as scr

def format_fasta(fasta, prefix = None, outdir = None):
    """
    Takes in the cellranger fasta files and add a prefix to the barcodes.
    This will generate a newly annotated fasta file that will be found in dandelion/data, unless otherwise specified in outdir

    Parameters
    ----------
    fasta
        fasta file
    prefix
        prefix to add
    outdir
        out location
    Returns
    -------
        new fasta file with new headers containing prefix
    """
    fh = open(fasta, 'r')
    seqs = {}
    for header, sequence in fasta_iterator(fh):
        if prefix is not None:
            newheader = prefix+'_'+str(header)
            seqs[newheader] = sequence
        else:
            seqs[header] = sequence
    fh.close()

    if os.path.isfile(fasta):
        basedir = os.path.dirname(fasta)
    elif os.path.isdir(fasta):
        basedir = os.path.dirname(fasta)
    else:
        basedir = os.getcwd()

    if outdir is None:
        out_dir = basedir+'/'+'dandelion/data/'
    else:
        if not outdir.endswith('/'):
            out_dir = outdir + '/'

    if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    out_fasta = out_dir + os.path.basename(fasta)

    fh1 = open(out_fasta, 'w')
    fh1.close()
    out = ''
    for l in seqs:
        out = '>'+l+'\n'+seqs[l]+'\n'
        Write_output(out, out_fasta)

    # format the barcode and contig_id in the corresponding annotation file too
    anno = basedir+'/'+os.path.basename(fasta).replace('.fasta', '_annotations.csv')
    data = pd.read_csv(anno, dtype = 'object')
    data['contig_id'] = [prefix+'_'+str(c) for c in data['contig_id']]
    data['barcode'] = [prefix+'_'+str(b).split('-')[0] for b in data['barcode']]
    out_anno = out_dir + os.path.basename(fasta).replace('.fasta', '_annotations.csv')
    data.to_csv(out_anno, index= False)    

def format_fastas(fastas, prefixes = None, outdir = None):
    """
    Takes in the cellranger fasta files and add a prefix to the barcodes.
    This will generate a newly annotated fasta file that will be found in dandelion/data, unless otherwise specified in outdir

    Parameters
    ----------
    fastas
        list or sequence of fasta files.
    prefixes
        list or sequence of prefixes to append to headers in each fasta file.
        if provided, a dictionary will be created with fasta files as keys and prefixes as values.
    outdir
        out location. Defaults to a subfolder called dandelion/data.
    Returns
    -------
        new fasta files with new headers containing prefixes
    """
    if type(fastas) is not list:
        fastas = [fastas]
    if prefixes is not None:
        if type(prefixes) is not list:
            prefixes = [prefixes]
        prefix_dict = dict(zip(fastas, prefixes))

    for fasta in tqdm(fastas, desc = 'Formating fasta(s) '):
        if prefixes is not None:
            format_fasta(fasta, prefix_dict[fasta], outdir)
        else:
            format_fasta(fasta, None, outdir)

def assign_isotype(fasta, fileformat = 'airr', org = 'human', blastdb = None, allele = False, parallel = True, dirs = None, verbose = False):
    """
    Takes in fasta files and annotate the parsed data table with the constant region.

    Parameters
    ----------
    fastas
        list or sequence of fasta files.
    format
        file format of parsed vdj results
    blastdb
        path to blast database.
    allele
        whether or not to return allele calls
    parallel
        whether or not to trigger parallelization
    dirs
        location for input (and also output)
    verbose
        whether or not to print the commands in terminal
    Returns
    -------
        new fasta files with new headers containing prefixes
    """
    def _run_blastn(fasta, blastdb, dirs, fileformat, org, verbose):

        env = os.environ.copy()
        if blastdb is None:
            try:
                bdb = env['BLASTDB']
            except:
                raise OSError('Environmental variable BLASTDB must be set. Otherwise, please provide path to blast database')
            bdb = bdb+org+'/'+org+'_BCR_C.fasta'
        else:
            env['BLASTDB'] = blastdb
            bdb = blastdb

        cmd = ['blastn',
                '-db', bdb,
                '-evalue', '0.001',
                '-max_target_seqs', '1',
                '-outfmt', '5',
                '-query', fasta]

        if dirs is None:
            blast_out = "{}/tmp/{}.xml".format(os.path.dirname(fasta), os.path.basename(fasta).split('.fasta')[0]+fileformat)
        else:
            blast_out = "{}/{}.xml".format(dirs, os.path.basename(fasta).split('.fasta')[0]+fileformat)
        if verbose:
            print('Running command: %s\n' % (' '.join(cmd)))
        with open(blast_out, 'w') as out:
            run(cmd, stdout = out, env = env)


    def _parse_BLAST(fasta, dirs, fileformat):
        '''
        Parses BLAST output from output files and writes formatted output to BLAST
        output summary files
        '''

        def split_blast_file(filename):
            '''
            code adapted from http://stackoverflow.com/questions/19575702/pythonhow-to-split-file-into-chunks-by-the-occurrence-of-the-header-word
            '''
            token = '<Iteration>'
            chunks = []
            current_chunk = []

            with open(filename) as fh:
                for line in fh:
                    line = line.rstrip()

                    if line.startswith(token) and current_chunk:
                        chunks.append(current_chunk[:])
                        current_chunk = []
                    if not line.startswith("Total queries"):
                        current_chunk.append(line)

                chunks.append(current_chunk)
            return (chunks)

        def extract_blast_info(line):
            line = line.split()[0]
            info = line.split(">")[1]
            info = info.split("<")[0]
            return (info)

        if dirs is None:
            input_file = "{}/tmp/{}.xml".format(os.path.dirname(fasta), os.path.basename(fasta).split('.fasta')[0]+fileformat)
            output_file = "{}/tmp/{}.blastsummary.txt".format(os.path.dirname(fasta), os.path.basename(fasta).split('.fasta')[0]+fileformat)
        else:
            input_file = "{}/{}.xml".format(dirs, os.path.basename(fasta).split('.fasta')[0]+fileformat)
            output_file = "{}/{}.blastsummary.txt".format(dirs, os.path.basename(fasta).split('.fasta')[0]+fileformat)

        with open(output_file, 'w') as outfile:
            outfile.write("------------------\n##{}##\n------------------\n\n#BCR#\n\n".format(fasta))
            # Split result file into chunks corresponding to results for each query sequence.
            if os.path.isfile(input_file):
                blast_result_chunks = split_blast_file(input_file)
                for chunk in blast_result_chunks:
                    message = False
                    for line_x in chunk:
                        line_x= line_x.strip()
                        if line_x.startswith("<Iteration_query-def>"):
                            line = line_x.split(">")[1]
                            blast_query_name = line.split("<")[0]
                        elif line_x.startswith("<Hsp_evalue>"):
                            evalue = extract_blast_info(line_x)
                            evalue = format(float(evalue), '.0e')
                        elif line_x.startswith("<Hit_accession>"):
                            C_segment = extract_blast_info(line_x)
                            if "C-REGION" or "CH1" in C_segment:
                                C_segment = C_segment.split("_")[0]
                        elif line_x.startswith("<Hsp_bit-score>"):
                            bit_score = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_query-from>"):
                            q_start = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_query-to>"):
                            q_end = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_hit-from>"):
                            s_start = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_hit-to>"):
                            s_end = extract_blast_info(line_x)
                        elif line_x.startswith("<Iteration_query-len>"):
                            query_length = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_align-len>"):
                            align_length = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_gaps>"):
                            gaps = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_identity>"):
                            identity = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_qseq>"):
                            c_qseq = extract_blast_info(line_x)
                        elif line_x.startswith("<Iteration_message>No hits found"):
                            message = True
                            out_string = "##{blast_query_name}##\nNo C segment found\n\n".format(
                                                                blast_query_name=blast_query_name)
                            outfile.write(out_string)
                        # Create output string when reaching end of BLAST
                        # iteration result (marked by </Iteration>) and write
                        # to BLAST summary file
                        elif line_x.startswith("</Iteration>") and message is not True:
                            identity_pro = float(identity)/int(align_length)*100
                            identity_pro = format(identity_pro, '.2f')
                            mismatches = int(align_length) - int(identity)
                            #Account for reversed sequences
                            if int(s_start) > int(s_end):
                                blast_query_name = "reversed|" + blast_query_name
                                x, y = int(q_start), int(q_end)
                                q_start = int(query_length) - y + 1
                                q_end = int(query_length) - x + 1
                                s_start, s_end = s_end, s_start
                            intro_string = "##{}##\nC segment:\t{}\n\n".format(
                                            blast_query_name, C_segment)
                            header_string = ("Segment\tquery_id\tsubject_id\t% identity\talignment length\t"
                                            "mismatches\tgap opens\tgaps\tq start\tq end\ts start\ts end\t"
                                            "evalue\tbit score\n")
                            out_string = ("C\t{blast_query_name}\t{C_segment}\t{identity_pro}\t{align_length}\t{mismatches}\tNA\t{gaps}\t{q_start}\t{q_end}\t{s_start}\t{s_end}\t{evalue}\t{bit_score}\t{q_seq}\n\n").format(
                                            blast_query_name=blast_query_name,
                                            C_segment=C_segment, identity_pro=identity_pro, align_length=align_length,
                                            evalue=evalue, mismatches=mismatches, gaps=gaps, q_start=q_start,
                                            q_end=q_end, s_start=s_start, s_end=s_end, bit_score=bit_score, q_seq = c_qseq)
                            string_to_write = intro_string + header_string + out_string
                            outfile.write(string_to_write)


    def _get_C(fasta, dirs, fileformat, allele = False, parallel = True):

        def _get_C_call(fasta, contig_name, dirs, fileformat, allele = False):
            if dirs is None:
                blast_summary_file = "{}/tmp/{}.blastsummary.txt".format(os.path.dirname(fasta), os.path.basename(fasta).split('.fasta')[0]+fileformat)
            else:
                blast_summary_file = "{}/{}.blastsummary.txt".format(dirs, os.path.basename(fasta).split('.fasta')[0]+fileformat)

            C_seq, C_gene, C_ident, C_eval, C_bitscore, C_qstart, C_qend = None, None, None, None, None, None, None
            with open(blast_summary_file, 'r') as input:
                for line in input:
                    if line.startswith("C\t{contig_name}".format(
                        contig_name=contig_name)) or line.startswith("C\treversed|{contig_name}".format(contig_name=contig_name)):
                        C_gene = line.split("\t")[2]
                        C_ident = line.split("\t")[3]
                        C_seq = line.split("\t")[14]
                        C_eval = line.split("\t")[12]
                        C_bitscore = line.split("\t")[13]
                        C_qstart = line.split("\t")[8]
                        C_qend = line.split("\t")[9]

                        if "_CH1" or "_C-REGION" in C_gene:
                            C_gene = C_gene.split("_")[0]
            if not allele:
                try:
                    C_gene = C_gene.split('*')[0]
                except:
                    pass
            
            C_call, C_identity, C_sequence, C_support, C_score, C_start, C_end, = {}, {}, {}, {}, {}, {}, {}
            C_call[contig_name] = C_gene
            C_identity[contig_name] = C_ident
            C_sequence[contig_name] = C_seq
            C_support[contig_name] = C_eval
            C_score[contig_name] = C_bitscore
            C_start[contig_name] = C_qstart
            C_end[contig_name] = C_qend

            return(C_sequence, C_call, C_identity, C_support, C_score, C_start, C_end)

        fh = open(fasta, 'r')
        contigs = []
        for header, sequence in fasta_iterator(fh):
            contigs.append(header)
        fh.close()

        if parallel:
            num_cores = multiprocessing.cpu_count()
            results = ()
            results = Parallel(n_jobs=num_cores)(delayed(_get_C_call)(fasta, c, dirs, fileformat, allele) for c in tqdm(contigs, desc = 'Retrieving contant region calls, parallelizing with ' + str(num_cores) + ' cpus '))                                    
            # transform list of dicts to dict
            seq, call, ident, support, score, start, end = {}, {}, {}, {}, {}, {}, {}
            for r in range(0, len(results)):                
                _seq, _call, _ident, _support, _score, _start, _end = results[r]
                seq.update(_seq)
                call.update(_call)
                ident.update(_ident)
                support.update(_support)
                score.update(_score)
                start.update(_start)
                end.update(_end)
        else:
            seq, call, ident, support, score, start, end = {}, {}, {}, {}, {}, {}, {}
            for c in tqdm(contigs, desc = 'Retrieving contant region calls '):
                seq[c], call[c], ident[c], support[c], score[c], start[c], end[c] = _get_C_call(fasta, c, dirs, fileformat, allele)[c]
        return(seq, call, ident, support, score, start, end)

    def _transfer_c(data, c_dict, colname):
        _data = load_data(data)
        if colname not in _data.columns:
            _data = _data.merge(pd.DataFrame.from_dict(c_dict, orient = 'index', columns = [colname]), left_index = True, right_index = True)
        else:
            _data[colname] = pd.Series(c_dict)
        return(_data)

    def _add_cell(data):
        _data = load_data(data)
        _data['cell_id'] = [c.split('_contig')[0].split('-')[0] for c in _data['sequence_id']]
        return(_data)

    format_dict = {'changeo':'_igblast_db-pass', 'airr':'_igblast_gap'}

    # running blast using blast
    _run_blastn(fasta, blastdb, dirs, format_dict[fileformat], org, verbose)
    # parsing output into a summary.txt file
    _parse_BLAST(fasta, dirs, format_dict[fileformat])
    # Add the c_calls to the data file
    c_seq, c_call, c_ident, c_supp, c_scr, c_st, c_en = {}, {}, {}, {}, {}, {}, {}
    c_seq, c_call, c_ident, c_supp, c_scr, c_st, c_en = _get_C(fasta, dirs, format_dict[fileformat], allele, parallel)
    if dirs is None:
        _file = "{}/{}.tsv".format(os.path.dirname(fasta), os.path.basename(fasta).split('.fasta')[0]+format_dict[fileformat])
    else:
        _file = "{}/{}.tsv".format(dirs, os.path.basename(fasta).split('.fasta')[0]+ format_dict[fileformat])
    dat = _transfer_c(_file, c_call, 'c_call')
    dat = _transfer_c(dat, c_seq, 'c_sequence_alignment')
    dat = _transfer_c(dat, c_st, 'c_sequence_start')
    dat = _transfer_c(dat, c_en, 'c_sequence_end')
    dat = _transfer_c(dat, c_scr, 'c_score')
    dat = _transfer_c(dat, c_ident, 'c_identity')
    dat = _transfer_c(dat, c_supp, 'c_support')
    dat = _add_cell(dat)
    dat.to_csv(_file, sep = '\t', index=False)

def reannotate_genes(data, igblast_db = None, germline = None, org ='human', loci = 'ig', fileformat = 'airr', dirs = None, filtered = False, extended = False, verbose = False, *args):
    """
    reannotate genes with igblastn and parses to data output

    Parameters
    ----------
    data
        list or sequence of fasta file locations, or folder name containing fasta files. if provided as a single string, it will first be converted to a list; this allows for the function to be run on single/multiple samples.
    dirs
        path to input files. will also determine folder structure for outout
    filtered
        whether or not the file being worked on is filtered or all (based on cellranger).
    verbose
        whether or not to print the command used in the terminal.
    fileformat: str (Default: 'airr')
        format of the output data. Default is 'airr'. changeo' will trigger legacy changeo format.
    *args
        passed to ``assigngenes_igblast`` and ``makedb_igblast``.
    """
    if type(data) is not list:
        data = [data]
    if dirs is None:
        path = 'dandelion/data/'
    else:
        if not dirs.endswith('/'):
            path = dirs + '/'
        else:
            path = dirs

    for s in tqdm(data, desc = 'Assigning genes '):
        if os.path.isfile(str(s)):
            filePath = s
        else:
            if filtered:
                filePath = s+'/'+path+'filtered_contig.fasta'
            else:
                filePath = s+'/'+path+'all_contig.fasta'
        assigngenes_igblast(filePath, igblast_db=igblast_db, org = org, loci=loci, fileformat = fileformat, outputfolder = dirs, verbose = verbose, *args)
        if fileformat == 'airr':
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
            insertGaps("{}/{}".format(os.path.dirname(filePath), os.path.basename(filePath).replace('.fasta', '_igblast.tsv')), [gml])
            map_cellranger("{}/{}".format(os.path.dirname(filePath), os.path.basename(filePath).replace('.fasta', '_igblast_gap.tsv')), extended = extended)
            tmpFolder = "{}/tmp".format(os.path.dirname(filePath))
            if not os.path.exists(tmpFolder):
                os.makedirs(tmpFolder)
            os.replace("{}/{}".format(os.path.dirname(filePath),os.path.basename(filePath).replace('.fasta', '_igblast.tsv')), "{}/{}".format(tmpFolder,os.path.basename(filePath).replace('.fasta', '_igblast.tsv')))            
        elif fileformat == 'changeo':            
            makedb_igblast(filePath, org = org, germline = germline, extended = extended, verbose = verbose)

def map_cellranger(data, extended = False):
    dat = load_data(data)
    cellranger_data = "{}/{}".format(os.path.dirname(data), os.path.basename(data).replace('_igblast_gap.tsv', '_annotations.csv'))
    cr_data = pd.read_csv(cellranger_data, dtype = 'object')
    cell_id = dict(zip(cr_data['contig_id'], cr_data['barcode']))
    v_call = dict(zip(cr_data['contig_id'], cr_data['v_gene']))
    d_call = dict(zip(cr_data['contig_id'], cr_data['d_gene']))
    j_call = dict(zip(cr_data['contig_id'], cr_data['j_gene']))
    c_call = dict(zip(cr_data['contig_id'], cr_data['c_gene']))
    junction = dict(zip(cr_data['contig_id'], cr_data['cdr3_nt']))
    junction_aa = dict(zip(cr_data['contig_id'], cr_data['cdr3']))
    conscount = dict(zip(cr_data['contig_id'], cr_data['reads']))
    umicount = dict(zip(cr_data['contig_id'], cr_data['umis']))

    if not extended:
        dat['cell_id'] = pd.Series(cell_id)
        dat['c_call'] = pd.Series(c_call)
        dat['consensus_count'] = pd.Series(conscount)
        dat['umi_count'] = pd.Series(umicount)
    else:
        dat['cell_id'] = pd.Series(cell_id)
        dat['c_call'] = pd.Series(c_call)
        dat['consensus_count'] = pd.Series(conscount)
        dat['umi_count'] = pd.Series(umicount)
        dat['v_call_10x'] = pd.Series(v_call)
        dat['d_call_10x'] = pd.Series(d_call)
        dat['j_call_10x'] = pd.Series(j_call)
        dat['junction_10x'] = pd.Series(junction)
        dat['junction_10x_aa'] = pd.Series(junction_aa)
    dat.to_csv(data, sep = '\t', index = False, na_rep='')

def reassign_alleles(data, out_folder, germline = None, fileformat = 'airr', plot = True, figsize = (4,3), dirs = None, sample_dict = None, filtered = False, out_filename = None, verbose = False):
    """
    Correct allele calls based on a personalized genotype.
    Description
    ----------
    reassignAlleles uses a subject-specific genotype to correct correct preliminary allele assignments of a set of sequences derived from a single subject.

    Parameters
    ----------
    data
        list or sequence of folders/data file locations. if provided as a single string, it will first be converted to a list; this allows for the function to be run on single/multiple samples.
    out_folder
        name of folder for concatenated data file and genotyped files.
    germline
        path to germline. If None, defaults to path set as environmental variable.
    fileformat
        format of data. only invoked if all other default options are used.
    dirs
        path to input files. will also determine folder structure for outout.
    sample_dict
        dictionary for creating a sample column
    filtered
        whether or not the file being worked on is filtered or all (based on cellranger).
    out_filename
        if provided, will save to this filename.
    verbose
        whether or not to print the command used in the terminal.
    *args
        passed to ``tigger_genotype``
    """
    def _return_IGKV_IGLV(results, locus = 'IGH'):
        res = results.copy()
        for i in tqdm(res.index, desc = '   Returning light chain V calls'):
            if ~(res.loc[i]['locus'] == locus):
                res.loc[i]['v_call_genotyped'] = res.loc[i]['v_call']
        return(res)

    if type(data) is not list:
        data = [data]
    if dirs is None:
        path = 'dandelion/data/'
    else:
        if not dirs.endswith('/'):
            path = dirs + '/'
        else:
            path = dirs

    if out_filename is not None:
        if not out_filename.endswith('.tsv'):
            raise OSError('Please provide a file name that ends with .tsv')

    informat_dict = {'changeo':'_igblast_db-pass.tsv', 'airr':'_igblast_gap.tsv'}
    fileformat_dict = {'changeo':'_igblast_db-pass_genotyped.tsv', 'airr':'_igblast_gap_genotyped.tsv'}
    inferred_fileformat_dict = {'changeo':'_igblast_db-pass_inferredGenotype.txt', 'airr':'_igblast_gap_inferredGenotype.txt'}
    data_list = []
    for s in tqdm(data, desc = 'Processing data file(s) '):
        if os.path.isfile(str(s)):
            filePath = s
        else:
            if filtered:
                filePath = s+'/'+path+'filtered_contig'+informat_dict[fileformat]
            else:
                filePath = s+'/'+path+'all_contig'+informat_dict[fileformat]
        dat = load_data(filePath)

        if sample_dict is not None:
            dat['sample_id'] = sample_dict[s]
        else:
            dat['sample_id'] = str(s)
        data_list.append(dat)

    # concatenate
    if len(data_list) > 1:
        print('Concatenating objects')
        dat_ = pd.concat(data_list, sort=False)
    else:
        dat_ = data_list[0]

    # dat_.fillna('', inplace=True)

    # write out this file for tigger
    if not out_folder.endswith('/'):
        outDir = out_folder + '/' + path
    else:
        outDir = out_folder + path

    if not os.path.exists(outDir):
        os.makedirs(outDir)
    if out_filename is None:
        if filtered:
            print('   Writing out concatenated object')
            # dat_.to_csv(outDir+'filtered_contig'+informat_dict[fileformat], index = False, sep = '\t', na_rep='')
            dat_h = dat_[dat_['locus'] == 'IGH']
            dat_h.to_csv(outDir+'filtered_contig_heavy'+informat_dict[fileformat], index = False, sep = '\t', na_rep='')
            tigger_genotype(outDir+'filtered_contig_heavy'+informat_dict[fileformat], germline = germline, fileformat = fileformat, verbose = verbose)
        else:
            print('   Writing out concatenated object')
            # dat_.to_csv(outDir+'all_contig'+informat_dict[fileformat], index = False, sep = '\t', na_rep='')
            dat_h = dat_[dat_['locus'] == 'IGH']
            dat_h.to_csv(outDir+'all_contig_heavy'+informat_dict[fileformat], index = False, sep = '\t', na_rep='')
            tigger_genotype(outDir+'all_contig_heavy'+informat_dict[fileformat], germline = germline, fileformat = fileformat, verbose = verbose)
    else:
        print('   Writing out concatenated object')
        # dat_.to_csv(out_filename, index = False, sep = '\t', na_rep='')
        dat_h = dat_[dat_['locus'] == 'IGH']
        dat_h.to_csv(outDir+'heavy_'+out_filename, index = False, sep = '\t', na_rep='')
        tigger_genotype(outDir+'heavy_'+out_filename, germline = germline, fileformat = fileformat, verbose = verbose)

    # and now to add it back to the original folders
    sleep(0.5)
    if out_filename is None:
        if filtered:
            out_h = load_data(outDir+'filtered_contig_heavy'+fileformat_dict[fileformat])
            # out = pd.read_csv(outDir+'filtered_contig'+fileformat_dict[fileformat], sep = '\t', dtype = 'object')
            dat_['v_call_genotyped'] = pd.Series(out_h['v_call_genotyped'])
            dat_ = _return_IGKV_IGLV(dat_)
            print('   Saving corrected genotyped object')
            sleep(0.5)
            dat_.to_csv(outDir+'filtered_contig'+fileformat_dict[fileformat], index = False, sep = '\t')
        else:
            out_h = load_data(outDir+'all_contig_heavy'+fileformat_dict[fileformat])
            # out = pd.read_csv(outDir+'all_contig'+fileformat_dict[fileformat], sep = '\t', dtype = 'object')
            dat_['v_call_genotyped'] = pd.Series(out_h['v_call_genotyped'])
            dat_ = _return_IGKV_IGLV(dat_)
            print('   Saving corrected genotyped object')
            sleep(0.5)
            dat_.to_csv(outDir+'all_contig'+fileformat_dict[fileformat], index = False, sep = '\t')
    else:
        out_h = load_data(outDir+'heavy_'+out_filename.replace('.tsv', '_genotyped.tsv'))
        # out = pd.read_csv(outDir+out_filename.replace('.tsv', '_genotyped.tsv'), sep = '\t', dtype = 'object')
        dat_['v_call_genotyped'] = pd.Series(out_h['v_call_genotyped'])
        dat_ = _return_IGKV_IGLV(dat_)
        print('   Saving corrected genotyped object')
        sleep(0.5)
        dat_.to_csv(out_filename.replace('.tsv', '_genotyped.tsv'), index = False, sep = '\t')

    for s in tqdm(data, desc = 'Writing out to individual folders '):
        if sample_dict is not None:
            out_ = dat_[dat_['sample_id'] == sample_dict[s]]
        else:
            out_ = dat_[dat_['sample_id'] == s]
        if os.path.isfile(str(s)):
            out_.to_csv(s.replace('.tsv', '_genotyped.tsv'), index = False, sep = '\t')
        else:
            if filtered:
                filePath = s+'/'+path+'filtered_contig'+fileformat_dict[fileformat]
            else:
                filePath = s+'/'+path+'all_contig'+fileformat_dict[fileformat]
            out_.to_csv(filePath, index = False, sep = '\t')
    if plot:
        print('Returning summary plot')
        if out_filename is None:
            if filtered:
                inferred_genotype = outDir+'filtered_contig_heavy'+inferred_fileformat_dict[fileformat]
            else:
                inferred_genotype = outDir+'all_contig_heavy'+inferred_fileformat_dict[fileformat]
        else:
            inferred_genotype = outDir+'heavy_'+out_filename.replace('.tsv', '_inferredGenotype.txt')
        inf_geno = pd.read_csv(inferred_genotype, sep = '\t', dtype = 'object')

        s2 = set(inf_geno['gene'])
        results = []
        for samp in list(set(out_h['sample_id'])):
            res = out_h[(out_h['sample_id']==samp)]
            V_ = [re.sub('[*][0-9][0-9]', '', v) for v in res['v_call']]
            V_g = [re.sub('[*][0-9][0-9]', '', v) for v in res['v_call_genotyped']]
            s1 = set(list(','.join([','.join(list(set(v.split(',')))) for v in V_]).split(',')))
            setdiff = s1 - s2
            ambiguous = (["," in i for i in V_].count(True)/len(V_)*100, ["," in i for i in V_g].count(True)/len(V_g)*100)
            not_in_genotype=([i in setdiff for i in V_].count(True)/len(V_)*100, [i in setdiff for i in V_g].count(True)/len(V_g)*100)
            stats = pd.DataFrame([ambiguous,not_in_genotype], columns = ['ambiguous', 'not_in_genotype'], index = ['before', 'after']).T
            stats.index.set_names(['vgroup'], inplace = True)
            stats.reset_index(drop = False, inplace = True)
            stats['sample_id'] = samp
            # stats['donor'] = str(out_folder)
            results.append(stats)
        results = pd.concat(results)
        ambiguous_table = results[results['vgroup'] == 'ambiguous']
        not_in_genotype_table = results[results['vgroup'] == 'not_in_genotype']
        ambiguous_table.reset_index(inplace = True, drop = True)
        not_in_genotype_table.reset_index(inplace = True, drop = True)
        # melting the dataframe
        ambiguous_table_before = ambiguous_table.drop('after', axis = 1)
        ambiguous_table_before.rename(columns={"before": "var"}, inplace = True)
        ambiguous_table_before['var_group'] = 'before'
        ambiguous_table_after = ambiguous_table.drop('before', axis = 1)
        ambiguous_table_after.rename(columns={"after": "var"}, inplace = True)
        ambiguous_table_after['var_group'] = 'after'
        ambiguous_table = pd.concat([ambiguous_table_before, ambiguous_table_after])
        not_in_genotype_table_before = not_in_genotype_table.drop('after', axis = 1)
        not_in_genotype_table_before.rename(columns={"before": "var"}, inplace = True)
        not_in_genotype_table_before['var_group'] = 'before'
        not_in_genotype_table_after = not_in_genotype_table.drop('before', axis = 1)
        not_in_genotype_table_after.rename(columns={"after": "var"}, inplace = True)
        not_in_genotype_table_after['var_group'] = 'after'
        not_in_genotype_table = pd.concat([not_in_genotype_table_before, not_in_genotype_table_after])
        ambiguous_table['var_group'] = ambiguous_table['var_group'].astype('category')
        not_in_genotype_table['var_group'] = not_in_genotype_table['var_group'].astype('category')
        ambiguous_table['var_group'].cat.reorder_categories(['before', 'after'], inplace = True)
        not_in_genotype_table['var_group'].cat.reorder_categories(['before', 'after'], inplace = True)

        options.figure_size = figsize
        final_table = pd.concat([ambiguous_table, not_in_genotype_table])
        p = (ggplot(final_table, aes('sample_id', y = 'var', fill='var_group'))
            + coord_flip()
            + theme_classic()
            + xlab("sample_id")
            + ylab("% allele calls")
            + ggtitle("Genotype reassignment with TIgGER")
            + geom_bar(stat="identity")
            + facet_grid('~'+str('vgroup'), scales="free_y")
            + scale_fill_manual(values=('#86bcb6', '#F28e2b'))
            + theme(legend_title = element_blank()))
        return(p)

def create_germlines(self, germline = None, org = 'human', seq_field='sequence_alignment', v_field='v_call', d_field='d_call', j_field='j_call', clone_field='clone_id', germ_types='dmask', fileformat='airr'):
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

    def _parseChangeO(record):
        """
        Parses a dictionary to a Receptor object

        Arguments:
          record : dict with fields and values in the Change-O format

        Returns:
          changeo.Receptor.Receptor : parsed Receptor object.
        """
        # Parse fields
        result = {}
        for k, v in record.items():
            k = ChangeoSchema.toReceptor(k)
            result[k] = v

        return Receptor(result)

    def _parseAIRR(record):
        """
        Parses a dictionary of AIRR records to a Receptor object

        Arguments:
          record : dict with fields and values in the AIRR format.

        Returns:
          changeo.Receptor.Receptor : parsed Receptor object.
        """
        # Parse fields
        result = {}
        for k, v in record.items():
            # Rename fields
            k = AIRRSchema.toReceptor(k)
            # Convert start positions to 0-based
            # if k in ReceptorData.start_fields and v is not None and v != '':
            #     v = str(int(v) + 1)
            # Assign new field
            result[k] = v

        for end, (start, length) in ReceptorData.end_fields.items():
            if end in result and result[end] is not None:
                try:
                    result[length] = int(result[end]) - int(result[start]) + 1
                except:
                    pass

        return Receptor(result)

    def _create_germlines_object(self, references, seq_field, v_field, d_field, j_field, clone_field, germ_types, fileformat):
        """
        Write germline sequences to tab-delimited database file

        Arguments:
        self : dandelion_class object
        references : folders and/or files containing germline repertoire data in FASTA format.
        seq_field : field in which to look for sequence.
        v_field : field in which to look for V call.
        d_field : field in which to look for D call.
        j_field : field in which to look for J call.
        cloned : if True build germlines by clone, otherwise build individual germlines.
        clone_field : field containing clone identifiers; ignored if cloned=False.
        germ_types : list of germline sequence types to be output from the set of 'full', 'dmask', 'vonly', 'regions'
        fileformat : input and output format.

        Returns:
        """
        # Define format operators
        try:
            reader, writer, schema = getFormatOperators(fileformat)
        except:
            raise ValueError('Invalid format %s' % fileformat)

        # Define output germline fields
        germline_fields = OrderedDict()
        seq_type = seq_field.split('_')[-1]
        if 'full' in germ_types:
            germline_fields['full'] = 'germline_' + seq_type
        if 'dmask' in germ_types:
            germline_fields['dmask'] = 'germline_' + seq_type + '_d_mask'
        if 'vonly' in germ_types:
            germline_fields['vonly'] = 'germline_' + seq_type + '_v_region'
        if 'regions' in germ_types:
            germline_fields['regions'] = 'germline_regions'

        if type(references) is not list:
            ref = [references]
        else:
            ref = ref
        reference_dict = readGermlines(ref)
        # Check for IMGT-gaps in germlines
        if all('...' not in x for x in reference_dict.values()):
            warnings.warn(UserWarning('Germline reference sequences do not appear to contain IMGT-numbering spacers. Results may be incorrect.'))

        required = ['v_germ_start_imgt', 'd_germ_start', 'j_germ_start', 'np1_length', 'np2_length']

        if self.__class__ == Dandelion:
            if isinstance(self.data, pd.DataFrame):
                # Check for required columns
                try:
                    checkFields(required, self.data.columns, schema=schema)
                except LookupError as e:
                    print(e)

                # Count input
                total_count = len(self.data)

                # Check for existence of fields
                for f in [v_field, d_field, j_field, seq_field]:
                    if f not in self.data.columns:
                        raise NameError('%s field does not exist in input database file.' % f)
                # Translate to Receptor attribute names
                v_field = schema.toReceptor(v_field)
                d_field = schema.toReceptor(d_field)
                j_field = schema.toReceptor(j_field)
                seq_field = schema.toReceptor(seq_field)
                clone_field = schema.toReceptor(clone_field)

                # Define Receptor iterator
                receptor_iter = ((self.data.loc[x, ].sequence_id, self.data.loc[x, ]) for x in self.data.index)
            else:
                raise LookupError('Please initialise the Dandelion object with a dataframe in data slot.')
        elif self.__class__ == pd.DataFrame:
            try:
                checkFields(required, self.columns, schema=schema)
            except LookupError as e:
                print(e)

            # Count input
            total_count = len(self)
            # Check for existence of fields
            for f in [v_field, d_field, j_field, seq_field]:
                if f not in self.columns:
                    raise NameError('%s field does not exist in input database file.' % f)
            # Translate to Receptor attribute names
            v_field = schema.toReceptor(v_field)
            d_field = schema.toReceptor(d_field)
            j_field = schema.toReceptor(j_field)
            seq_field = schema.toReceptor(seq_field)
            clone_field = schema.toReceptor(clone_field)
            # Define Receptor iterator
            receptor_iter = ((self.loc[x, ].sequence_id, self.loc[x, ]) for x in self.index)
        
        out = {}
        # Iterate over rows
        for key, records in tqdm(receptor_iter, desc = 'Building germline sequences '):
            # Define iteration variables
            # Build germline for records
            if fileformat == 'airr':
                germ_log, glines, genes = buildGermline(_parseAIRR(dict(records)), reference_dict, seq_field=seq_field, v_field=v_field, d_field=d_field, j_field=j_field)
            elif fileformat == 'changeo':
                germ_log, glines, genes = buildGermline(_parseChangeO(dict(records)), reference_dict, seq_field=seq_field, v_field=v_field, d_field=d_field, j_field=j_field)
            else:
                raise AttributeError('%s is not acceptable file format.' % fileformat)
            
            if glines is not None:
                # Add glines to Receptor record
                annotations = {}
                if 'full' in germ_types:
                    annotations[germline_fields['full']] = glines['full']
                if 'dmask' in germ_types:
                    annotations[germline_fields['dmask']] = glines['dmask']
                if 'vonly' in germ_types:
                    annotations[germline_fields['vonly']] = glines['vonly']
                if 'regions' in germ_types:
                    annotations[germline_fields['regions']] = glines['regions']
                out.update({key:annotations})
        germline_df = pd.DataFrame.from_dict(out, orient = 'index')

        if self.__class__ == Dandelion:
            self.data = load_data(self.data)
            for x in germline_df.columns:
                self.data[x] = pd.Series(germline_df[x])
        elif self.__class__ == pd.DataFrame:
            out = Dandelion(self)
            for x in germline_df.columns:
                out.data[x] = pd.Series(germline_df[x])
            return(out)
        
    def _create_germlines_file(file, references, seq_field, v_field, d_field, j_field, clone_field, germ_types, fileformat):
        """
        Write germline sequences to tab-delimited database file

        Arguments:
        file : airr/changeo tsv file
        references : folders and/or files containing germline repertoire data in FASTA format.
        seq_field : field in which to look for sequence.
        v_field : field in which to look for V call.
        d_field : field in which to look for D call.
        j_field : field in which to look for J call.
        cloned : if True build germlines by clone, otherwise build individual germlines.
        clone_field : field containing clone identifiers; ignored if cloned=False.
        germ_types : list of germline sequence types to be output from the set of 'full', 'dmask', 'vonly', 'regions'
        fileformat : input and output format.

        Returns:
        """
        # Define format operators
        try:
            reader, writer, schema = getFormatOperators(fileformat)
        except:
            raise ValueError('Invalid format %s' % fileformat)

        # Define output germline fields
        germline_fields = OrderedDict()
        seq_type = seq_field.split('_')[-1]
        if 'full' in germ_types:
            germline_fields['full'] = 'germline_' + seq_type
        if 'dmask' in germ_types:
            germline_fields['dmask'] = 'germline_' + seq_type + '_d_mask'
        if 'vonly' in germ_types:
            germline_fields['vonly'] = 'germline_' + seq_type + '_v_region'
        if 'regions' in germ_types:
            germline_fields['regions'] = 'germline_regions'

        if type(references) is not list:
            ref = [references]
        else:
            ref = ref
        reference_dict = readGermlines(ref)
        # Check for IMGT-gaps in germlines
        if all('...' not in x for x in reference_dict.values()):
            warnings.warn(UserWarning('Germline reference sequences do not appear to contain IMGT-numbering spacers. Results may be incorrect.'))

        required = ['v_germ_start_imgt', 'd_germ_start', 'j_germ_start', 'np1_length', 'np2_length']

        # Get repertoire and open Db reader
        db_handle = open(file, 'rt')
        db_iter = reader(db_handle)
        # Check for required columns
        try:
            checkFields(required, db_iter.fields, schema=schema)
        except LookupError as e:
            print(e)
        # Count input
        total_count = countDbFile(file)
        # Check for existence of fields
        for f in [v_field, d_field, j_field, seq_field]:
            if f not in db_iter.fields:
                raise NameError('%s field does not exist in input database file.' % f)
        # Translate to Receptor attribute names
        v_field = schema.toReceptor(v_field)
        d_field = schema.toReceptor(d_field)
        j_field = schema.toReceptor(j_field)
        seq_field = schema.toReceptor(seq_field)
        clone_field = schema.toReceptor(clone_field)
        # Define Receptor iterator
        receptor_iter = ((x.sequence_id, [x]) for x in db_iter)
        
        out = {}
        # Iterate over rows
        for key, records in tqdm(receptor_iter, desc = 'Building germline sequences '):
            # Define iteration variables
            # Build germline for records
            # if not isinstance(self.data, pd.DataFrame):
            records = list(records)
            germ_log, glines, genes = buildGermline(records[0], reference_dict, seq_field=seq_field, v_field=v_field, d_field=d_field, j_field=j_field)
            if glines is not None:
                # Add glines to Receptor record
                annotations = {}
                if 'full' in germ_types:
                    annotations[germline_fields['full']] = glines['full']
                if 'dmask' in germ_types:
                    annotations[germline_fields['dmask']] = glines['dmask']
                if 'vonly' in germ_types:
                    annotations[germline_fields['vonly']] = glines['vonly']
                if 'regions' in germ_types:
                    annotations[germline_fields['regions']] = glines['regions']
                out.update({key:annotations})
        germline_df = pd.DataFrame.from_dict(out, orient = 'index')

        out = Dandelion(data = file)
        for x in germline_df.columns:
            out.data[x] = pd.Series(germline_df[x])

        if os.path.isfile(str(file)):
            out.data.to_csv("{}/{}_germline_{}.tsv".format(os.path.dirname(file), os.path.basename(file).split('.tsv')[0], germ_types), sep = '\t', index = False)
        
        return(out)

    if self.__class__ == Dandelion or self.__class__ == pd.DataFrame:
        _create_germlines_object(self, gml, seq_field, v_field, d_field, j_field, clone_field, germ_types, fileformat)
    else:
        return(_create_germlines_file(self, gml, seq_field, v_field, d_field, j_field, clone_field, germ_types, fileformat))

def run_scanpy_qc(self, max_genes=2500, min_genes=200, mito_cutoff=0.05, pval_cutoff=0.1, min_counts=None, max_counts=None):
    """
    Parameters
    ----------
    adata : AnnData
        The (annotated) data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    max_genes : int (default: 2500)
        Maximum number of genes expressed required for a cell to pass filtering.
    min_genes : int (default: 200)
        Minimum number of genes expressed  required for a cell to pass filtering.
    mito_cutoff : float (default: 0.05)
        Maximum percentage mitochondrial content allowed for a cell to pass filtering.
    pval_cutoff : float (default: 0.05)
        Maximum Benjamini-Hochberg corrected p value from doublet detection protocol allowed for a cell to pass filtering.
    min_counts : int, None (default: None)
        Minimum number of counts required for a cell to pass filtering.
    max_counts : int, None (default: None)
        Maximum number of counts required for a cell to pass filtering.
    Returns
    -------
    AnnData
        The (annotated) data matrix of shape n_obs × n_vars where obs now contain filtering information. Rows correspond to cells and columns to genes.

    """
    _adata = self.copy()
    # run scrublet
    scrub = scr.Scrublet(_adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    _adata.obs['scrublet_score'] = doublet_scores
    # overcluster prep. run basic scanpy pipeline
    sc.pp.filter_cells(_adata, min_genes = 0)
    mito_genes = _adata.var_names.str.startswith('MT-')
    _adata.obs['percent_mito'] = np.sum(_adata[:, mito_genes].X, axis = 1) / np.sum(_adata.X, axis = 1)
    _adata.obs['n_counts'] = _adata.X.sum(axis = 1)
    sc.pp.normalize_total(_adata)
    sc.pp.log1p(_adata)
    sc.pp.highly_variable_genes(_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
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
    # let's do a one-sided test. the Bertie write-up does not address this but it makes sense
    pvals = 1-scipy.stats.norm.cdf(_adata.obs['scrublet_cluster_score'], loc=med, scale=1.4826*mad)
    _adata.obs['bh_pval'] = bh(pvals)
    # threshold the p-values to get doublet calls.
    _adata.obs['is_doublet'] = _adata.obs['bh_pval'] < pval_cutoff
    _adata.obs['is_doublet'] = _adata.obs['is_doublet'].astype('category')
    _adata.obs['filter_rna'] = (pd.Series([min_genes < n > max_genes for n in _adata.obs['n_genes']], index = _adata.obs.index)) | \
        (_adata.obs['percent_mito'] >= mito_cutoff) | \
            (_adata.obs['is_doublet'] == True)

    # removing columns that probably don't need anymore
    _adata.obs = _adata.obs.drop(['leiden', 'leiden_R', 'scrublet_cluster_score'], axis = 1)
    self.obs = _adata.obs.copy()