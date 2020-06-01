#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 17:56:02
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-06-01 22:12:41

import sys
import os
import pandas as pd
from subprocess import run
from tqdm import tqdm
import multiprocessing
from joblib import Parallel, delayed
from time import sleep
from ..utilities._misc import *
from .ext._immcantationscripts import assigngenes_igblast, makedb_igblast, tigger_genotype, insertGaps

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
            bdb = env['BLASTDB']

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
                            out_string = ("C\t{blast_query_name}\t{C_segment}\t{identity_pro}\t{align_length}\t{mismatches}\tNA\t{gaps}\t{q_start}\t{q_end}\t{s_start}\t{s_end}\t{evalue}\t{bit_score}\n\n").format(
                                            blast_query_name=blast_query_name,
                                            C_segment=C_segment, identity_pro=identity_pro, align_length=align_length,
                                            evalue=evalue, mismatches=mismatches, gaps=gaps, q_start=q_start,
                                            q_end=q_end, s_start=s_start, s_end=s_end, bit_score=bit_score)
                            string_to_write = intro_string + header_string + out_string
                            outfile.write(string_to_write)


    def _get_C(fasta, dirs, fileformat, allele = False, parallel = True):

        def _get_C_call(fasta, contig_name, dirs, fileformat, allele = False):
            if dirs is None:
                blast_summary_file = "{}/tmp/{}.blastsummary.txt".format(os.path.dirname(fasta), os.path.basename(fasta).split('.fasta')[0]+fileformat)
            else:
                blast_summary_file = "{}/{}.blastsummary.txt".format(dirs, os.path.basename(fasta).split('.fasta')[0]+fileformat)

            C_gene = None
            with open(blast_summary_file, 'r') as input:
                for line in input:
                    if line.startswith("C\t{contig_name}".format(
                        contig_name=contig_name)) or line.startswith("C\treversed|{contig_name}".format(contig_name=contig_name)):
                        C_gene = line.split("\t")[2]

                        if "_CH1" or "_C-REGION" in C_gene:
                            C_gene = C_gene.split("_")[0]
            if not allele:
                try:
                    C_gene = C_gene.split('*')[0]
                except:
                    pass
            C_call = {}
            C_call[contig_name] = C_gene
            return(C_call)

        fh = open(fasta, 'r')
        contigs = []
        for header, sequence in fasta_iterator(fh):
            contigs.append(header)
        fh.close()

        if parallel:
            num_cores = multiprocessing.cpu_count()
            results = {}
            results = Parallel(n_jobs=num_cores)(delayed(_get_C_call)(fasta, c, dirs, fileformat, allele) for c in tqdm(contigs, desc = 'Retrieving contant region calls, parallelizing with ' + str(num_cores) + ' cpus '))
            # transform list of dicts to dict
            results = {k: v for x in results for k, v in x.items()}
        else:
            results = {}
            for c in tqdm(contigs, desc = 'Retrieving contant region calls '):
                results[c] = _get_C_call(fasta, c, dirs, fileformat, allele)[c]
        return(results)

    def _transfer_c_call(data, c_dict):
        _data = load_data(data)
        if 'c_call' not in _data.columns:
            _data = _data.merge(pd.DataFrame.from_dict(c_dict, orient = 'index', columns = ['c_call']), left_index = True, right_index = True)
        else:
            _data['c_call'] = pd.Series(c_dict)
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
    c_call = _get_C(fasta, dirs, format_dict[fileformat])
    if dirs is None:
        _file = "{}/{}.tsv".format(os.path.dirname(fasta), os.path.basename(fasta).split('.fasta')[0]+format_dict[fileformat])
    else:
        _file = "{}/{}.tsv".format(dirs, os.path.basename(fasta).split('.fasta')[0]+ format_dict[fileformat])
    dat = _transfer_c_call(_file, c_call)
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

def reassign_alleles(data, out_folder, fileformat = 'airr', dirs = None, sample_dict = None, filtered = False, out_filename = None, verbose = False, *args):
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
        for i in tqdm(results.index, desc = '   Returning light chain V calls'):
            if ~(results.iloc[i]['locus'] == locus):
                results.iloc[i]['v_call_genotyped'] = results.iloc[i]['v_call']
        return(results)

    if type(data) is not list:
        data = [data]
    if dirs is None:
        path = 'dandelion/data/'
    else:
        if not dirs.endswith('/'):
            path = dirs + '/'
        else:
            path = dirs

    informat_dict = {'changeo':'_igblast_db-pass.tsv', 'airr':'_igblast_gap.tsv'}
    fileformat_dict = {'changeo':'_igblast_db-pass_genotyped.tsv', 'airr':'_igblast_gap_genotyped.tsv'}
    data_list = []
    for s in tqdm(data, desc = 'Processing data file(s) '):
        if os.path.isfile(str(s)):
            filePath = s
        else:
            if filtered:
                filePath = s+'/'+path+'filtered_contig'+informat_dict[fileformat]
            else:
                filePath = s+'/'+path+'all_contig'+informat_dict[fileformat]
        dat = pd.read_csv(filePath, sep = '\t', dtype = 'object')

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

    dat_.fillna('', inplace=True)

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
            dat_.to_csv(outDir+'filtered_contig'+informat_dict[fileformat], index = False, sep = '\t', na_rep='')
            tigger_genotype(outDir+'filtered_contig'+informat_dict[fileformat], fileformat = fileformat, verbose = verbose, *args)
        else:
            print('   Writing out concatenated object')
            dat_.to_csv(outDir+'all_contig'+informat_dict[fileformat], index = False, sep = '\t', na_rep='')
            tigger_genotype(outDir+'all_contig'+informat_dict[fileformat], fileformat = fileformat, verbose = verbose, *args)
    else:
        print('   Writing out concatenated object')
        dat_.to_csv(out_filename, index = False, sep = '\t', na_rep='')
        tigger_genotype(out_filename, verbose = verbose, fileformat = fileformat,  *args)

    # and now to add it back to the original folders
    print('   Reading genotyped object')
    sleep(0.5)
    if out_filename is None:
        if filtered:
            out = pd.read_csv(outDir+'filtered_contig'+fileformat_dict[fileformat], sep = '\t', dtype = 'object')
            out = _return_IGKV_IGLV(out)
            print('   Saving corrected genotyped object')
            sleep(0.5)
            out.to_csv(outDir+'filtered_contig'+fileformat_dict[fileformat], index = False, sep = '\t')
        else:
            out = pd.read_csv(outDir+'all_contig'+fileformat_dict[fileformat], sep = '\t', dtype = 'object')
            out = _return_IGKV_IGLV(out)
            print('   Saving corrected genotyped object')
            sleep(0.5)
            out.to_csv(outDir+'all_contig'+fileformat_dict[fileformat], index = False, sep = '\t')
    else:
        out = pd.read_csv(out_filename.replace('.tsv', '_genotyped.tsv'), sep = '\t', dtype = 'object')
        out = _return_IGKV_IGLV(out)
        print('   Saving corrected genotyped object')
        sleep(0.5)
        out.to_csv(out_filename.replace('.tsv', '_genotyped.tsv'), index = False, sep = '\t')

    for s in tqdm(data, desc = 'Writing out to individual folders '):
        if sample_dict is not None:
            out_ = out[out['sample_id'] == sample_dict[s]]
        else:
            out_ = out[out['sample_id'] == s]
        if os.path.isfile(str(s)):
            out_.to_csv(s.replace('.tsv', '_genotyped.tsv'), index = False, sep = '\t')
        else:
            if filtered:
                filePath = s+'/'+path+'filtered_contig'+fileformat_dict[fileformat]
            else:
                filePath = s+'/'+path+'all_contig'+fileformat_dict[fileformat]
            out_.to_csv(filePath, index = False, sep = '\t')