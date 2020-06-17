#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 17:56:02
# @Last Modified by:   Kelvin
# @Last Modified time: 2020-06-17 17:21:15

import sys
import os
import pandas as pd
from subprocess import run
from tqdm import tqdm
import multiprocessing
from joblib import Parallel, delayed
from collections import OrderedDict
from time import sleep
from ..utilities._utilities import *
from .ext._preprocessing import assigngenes_igblast, makedb_igblast, tigger_genotype, insertGaps
from plotnine import ggplot, geom_bar, ggtitle, scale_fill_manual, coord_flip, options, element_blank, aes, xlab, ylab, facet_grid, theme_classic, theme
from changeo.Gene import buildGermline
from changeo.IO import countDbFile, getDbFields, getFormatOperators, readGermlines, checkFields
from changeo.Receptor import AIRRSchema, ChangeoSchema, Receptor, ReceptorData
import re
import scanpy as sc
import numpy as np
import scipy.stats
import scrublet as scr
from Bio import Align

def format_fasta(fasta, prefix = None, outdir = None):
    """
    Adds prefix to the headers/contig ids in cellranger fasta and annotation file.

    Parameters
    ----------
    fasta : str
        path to fasta file.
    prefix : str, optional
        prefix to append to the headers/contig ids.
    outdir : str, optional
        path to out put location. Default is None, which is 'dandelion/data'.
    Returns
    -------
        Formatted fasta file with new headers containing prefix
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
    Adds prefix to the headers/contig ids in cellranger fasta and annotation file.

    Parameters
    ----------
    fastas : list
        list or sequence of paths to fasta files.
    prefixes : list, optional
        list or sequence of prefixes to append to headers/contig ids in each fasta file.
    outdir : str, optional
        path to out put location. Default is None, which is 'dandelion/data'.
    Returns
    -------
        Formatted fasta file with new headers containing prefix
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

def assign_isotype(fasta, fileformat = 'airr', org = 'human', correct_c_call = True, correction_dict = None, plot = True, figsize=(4,4), blastdb = None, allele = False, parallel = True, dirs = None, verbose = False):
    """
    Annotate contigs with constant region call using blastn

    Parameters
    ----------
    fasta : str
        path to fasta file.
    fileformat : str
        format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.
    org : str
        organism of reference folder. Default is 'human'.
    correct_c_call : bool
        whether or not to adjust the c_calls after blast based on provided primers specified in `primer_dict` option. Default is True.
    correction_dict : dict[dict], optional
        a nested dictionary contain isotype/c_genes as keys and primer sequences as records to use for correcting annotated c_calls. Defaults to a curated dictionary for human sequences if left as none.
    plot : bool
        whether or not to plot reassignment summary metrics. Default is True.
    figsize : tuple[float, float]
        size of figure. Default is (4, 4).
    blastdb : str, optional
        path to blast database. Defaults to `$BLASTDB` environmental variable.
    allele : bool
        whether or not to return allele calls. Default is False.
    parallel : bool
        whether or not to use parallelization. Default is True.
    dirs : str, optional
        location of both input and output files. None defaults to dandelion/data folder.
    verbose : bool
        whether or not to print the blast command in terminal. Default is False.
    Returns
    -------
        V(D)J tsv files with constant genes annotated.
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
                        elif line_x.startswith("<Hsp_hseq>"):
                            c_hseq = extract_blast_info(line_x)
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
                            out_string = ("C\t{blast_query_name}\t{C_segment}\t{identity_pro}\t{align_length}\t{mismatches}\tNA\t{gaps}\t{q_start}\t{q_end}\t{s_start}\t{s_end}\t{evalue}\t{bit_score}\t{q_seq}\t{h_seq}\n\n").format(
                                            blast_query_name=blast_query_name,
                                            C_segment=C_segment, identity_pro=identity_pro, align_length=align_length,
                                            evalue=evalue, mismatches=mismatches, gaps=gaps, q_start=q_start,
                                            q_end=q_end, s_start=s_start, s_end=s_end, bit_score=bit_score, q_seq = c_qseq, h_seq = c_hseq)
                            string_to_write = intro_string + header_string + out_string
                            outfile.write(string_to_write)


    def _get_C(fasta, dirs, fileformat, allele = False, parallel = True):

        def _get_C_call(fasta, contig_name, dirs, fileformat, allele = False):
            if dirs is None:
                blast_summary_file = "{}/tmp/{}.blastsummary.txt".format(os.path.dirname(fasta), os.path.basename(fasta).split('.fasta')[0]+fileformat)
            else:
                blast_summary_file = "{}/{}.blastsummary.txt".format(dirs, os.path.basename(fasta).split('.fasta')[0]+fileformat)

            C_seq,C_germ, C_gene, C_ident, C_eval, C_bitscore, C_qstart, C_qend = None, None, None, None, None, None, None, None
            with open(blast_summary_file, 'r') as input:
                for line in input:
                    if line.startswith("C\t{contig_name}".format(
                        contig_name=contig_name)) or line.startswith("C\treversed|{contig_name}".format(contig_name=contig_name)):
                        C_gene = line.split("\t")[2]
                        C_ident = line.split("\t")[3]
                        C_seq = line.split("\t")[14]
                        C_germ = line.split("\t")[15]
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
            
            C_call, C_identity, C_sequence, C_germline, C_support, C_score, C_start, C_end, = {}, {}, {}, {}, {}, {}, {}, {}
            C_call[contig_name] = C_gene
            C_identity[contig_name] = C_ident
            C_sequence[contig_name] = C_seq
            C_germline[contig_name] = C_germ
            C_support[contig_name] = C_eval
            C_score[contig_name] = C_bitscore
            C_start[contig_name] = C_qstart
            C_end[contig_name] = C_qend

            return(C_sequence, C_germline, C_call, C_identity, C_support, C_score, C_start, C_end)

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
            seq, germ, call, ident, support, score, start, end = {}, {}, {}, {}, {}, {}, {}, {}
            for r in range(0, len(results)):                
                _seq, _germ, _call, _ident, _support, _score, _start, _end = results[r]
                seq.update(_seq)
                germ.update(_germ)
                call.update(_call)
                ident.update(_ident)
                support.update(_support)
                score.update(_score)
                start.update(_start)
                end.update(_end)
        else:
            seq, germ, call, ident, support, score, start, end = {}, {}, {}, {}, {}, {}, {}, {}
            for c in tqdm(contigs, desc = 'Retrieving contant region calls '):
                seq[c], germ[c], call[c], ident[c], support[c], score[c], start[c], end[c] = _get_C_call(fasta, c, dirs, fileformat, allele)[c]
        return(seq, germ, call, ident, support, score, start, end)

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

    aligner = Align.PairwiseAligner()

    def two_gene_correction(self, i, dictionary):    
        key1, key2 = dictionary.keys()
        seq = self.loc[i, 'c_sequence_alignment'].replace('-', '')
        alignments1 = aligner.align(dictionary[key1], seq)
        alignments2 = aligner.align(dictionary[key2], seq)
        score1 = alignments1.score
        score2 = alignments2.score
        if score1 == score2:
            self.loc[i, 'c_call'] = str(key1)+','+str(key2)
        if score1 > score2:
            self.loc[i, 'c_call'] = str(key1)
        if score1 < score2:
            self.loc[i, 'c_call'] = str(key2)

    def three_gene_correction(self, i, dictionary):
        key1, key2, key3 = dictionary.keys()
        seq = self.loc[i, 'c_sequence_alignment'].replace('-', '')
        alignments1 = aligner.align(dictionary[key1], seq)
        alignments2 = aligner.align(dictionary[key2], seq)
        alignments3 = aligner.align(dictionary[key3], seq)
        score1 = alignments1.score
        score2 = alignments2.score
        score3 = alignments3.score
        if score1 == score2 == score3:
            self.loc[i, 'c_call'] = str(key1)+','+str(key2)+','+str(key3)
        elif score1 > score2 and score1 > score3:
            self.loc[i, 'c_call'] = str(key1)
        elif score2 > score1 and score2 > score3:
            self.loc[i, 'c_call'] = str(key2)
        elif score3 > score1 and score3 > score2:
            self.loc[i, 'c_call'] = str(key3)
        elif score1 == score2 and score1 > score3:
            self.loc[i, 'c_call'] = str(key1)+','+str(key2)
        elif score1 > score2 and score1 == score3:
            self.loc[i, 'c_call'] = str(key1)+','+str(key3)
        elif score2 > score1 and score2 == score3:
            self.loc[i, 'c_call'] = str(key2)+','+str(key3)
                
    def four_gene_correction(self, i, dictionary):    
        key1, key2, key3, key4 = dictionary.keys()
        seq = self.loc[i, 'c_sequence_alignment'].replace('-', '')
        alignments1 = aligner.align(dictionary[key1], seq)
        alignments2 = aligner.align(dictionary[key2], seq)
        alignments3 = aligner.align(dictionary[key3], seq)
        alignments4 = aligner.align(dictionary[key4], seq)
        score1 = alignments1.score
        score2 = alignments2.score
        score3 = alignments3.score
        score4 = alignments4.score
        if score1 == score2 == score3 == score4:
            self.loc[i, 'c_call'] = str(key1)+','+str(key2)+','+str(key3)+','+str(key4)
        elif score1 > score2 and score1 > score3 and score1 > score4:
            self.loc[i, 'c_call'] = str(key1)
        elif score2 > score1 and score2 > score3 and score2 > score4:
            self.loc[i, 'c_call'] = str(key2)
        elif score3 > score1 and score3 > score2 and score3 > score4:
            self.loc[i, 'c_call'] = str(key3)
        elif score4 > score1 and score4 > score2 and score4 > score3:
            self.loc[i, 'c_call'] = str(key4)
        elif score1 == score2 and score1 > score3 and score1 > score4:
            self.loc[i, 'c_call'] = str(key1)+','+str(key2)
        elif score1 > score2 and score1 == score3 and score1 > score4:
            self.loc[i, 'c_call'] = str(key1)+','+str(key3)
        elif score1 > score2 and score1 > score3 and score1 == score4:
            self.loc[i, 'c_call'] = str(key1)+','+str(key4)
        elif score2 == score3 and score2 > score1 and score2 > score4:
            self.loc[i, 'c_call'] = str(key1)+','+str(key3)
        elif score2 == score4 and score2 > score1 and score2 > score3:
            self.loc[i, 'c_call'] = str(key2)+','+str(key4)
        elif score3 == score4 and score3 > score1 and score3 > score2:
            self.loc[i, 'c_call'] = str(key3)+','+str(key4)
        elif score1 == score2 == score3 and score1 > score4:
            self.loc[i, 'c_call'] = str(key1)+','+str(key2)+','+str(key3)
        elif score1 == score2 == score4 and score1 > score3:
            self.loc[i, 'c_call'] = str(key1)+','+str(key2)+','+str(key4)
        elif score1 == score3 == score4 and score1 > score2:
            self.loc[i, 'c_call'] = str(key1)+','+str(key3)+','+str(key4)
        elif score2 == score3 == score4 and score2 > score1:
            self.loc[i, 'c_call'] = str(key2)+','+str(key3)+','+str(key4)
            
    def _correct_c_call(data, primers_dict=None):
        dat = data.copy()
        if primers_dict is None:
            primer_dict = {
                'IGHG':{
                    'IGHG1':'GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCCTGGGC',
                    'IGHG2':'GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGC',
                    'IGHG3':'GCTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGGGCACAGCGGCCCTGGGC',
                    'IGHG4':'GCTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGGC'},
                'IGHA':{
                    'IGHA1':'GCATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCAGCACCCAGCCAGATGGGAACGTGGTCATCGCCTGC',
                    'IGHA2':'GCATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACAGCACCCCCCAAGATGGGAACGTGGTCGTCGCATGC'},
                'IGLC7':{
                    'IGLC':'GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAA',
                    'IGLC7':'GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCACCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCGTAA'},
                'IGLC3':{
                    'IGLC':'GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAA',
                    'IGLC3':'GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCACCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAA'},
                'IGLC6':{
                    'IGLC': 'TCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCA',
                    'IGLC6':'TCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGCCTGA'}}
        else:
            primer_dict = primers_dict
        
        for i in dat.index:
            if (dat.loc[i, 'c_call'] is not np.nan) & (dat.loc[i, 'c_call'] is not None):
                for k in primer_dict:
                    if k in dat.loc[i, 'c_call']:
                        if len(primer_dict[k]) == 2:
                            two_gene_correction(dat, i, primer_dict[k])
                        elif len(primer_dict[k]) == 3:
                            three_gene_correction(dat, i, primer_dict[k])                        
                        elif len(primer_dict[k]) == 4:
                            four_gene_correction(dat, i, primer_dict[k])
        return(dat)

    format_dict = {'changeo':'_igblast_db-pass', 'airr':'_igblast_gap'}

    # running blast using blast
    _run_blastn(fasta, blastdb, dirs, format_dict[fileformat], org, verbose)
    # parsing output into a summary.txt file
    _parse_BLAST(fasta, dirs, format_dict[fileformat])
    # Add the c_calls to the data file
    c_seq, c_germ, c_call, c_ident, c_supp, c_scr, c_st, c_en = {}, {}, {}, {}, {}, {}, {}, {}
    c_seq, c_germ, c_call, c_ident, c_supp, c_scr, c_st, c_en = _get_C(fasta, dirs, format_dict[fileformat], allele, parallel)
    if dirs is None:
        _file = "{}/{}.tsv".format(os.path.dirname(fasta), os.path.basename(fasta).split('.fasta')[0]+format_dict[fileformat])
    else:
        _file = "{}/{}.tsv".format(dirs, os.path.basename(fasta).split('.fasta')[0]+ format_dict[fileformat])
    dat_10x = load_data(_file)
    res_10x = pd.DataFrame(dat_10x['c_call'])
    dat = _transfer_c(_file, c_call, 'c_call')
    dat = _transfer_c(dat, c_seq, 'c_sequence_alignment')
    dat = _transfer_c(dat, c_germ, 'c_germline_alignment')
    dat = _transfer_c(dat, c_st, 'c_sequence_start')
    dat = _transfer_c(dat, c_en, 'c_sequence_end')
    dat = _transfer_c(dat, c_scr, 'c_score')
    dat = _transfer_c(dat, c_ident, 'c_identity')
    dat = _transfer_c(dat, c_supp, 'c_support')
    res_blast = pd.DataFrame(dat['c_call'])  

    res_10x_sum = pd.DataFrame(res_10x['c_call'].value_counts(normalize=True)*100)
    res_blast_sum = pd.DataFrame(res_blast['c_call'].value_counts(normalize=True)*100)
    res_10x_sum['group'] = '10X'
    res_blast_sum['group'] = 'blast'
    res_10x_sum.columns = ['counts', 'group']
    res_blast_sum.columns = ['counts', 'group']
    res_10x_sum.index = res_10x_sum.index.set_names(['c_call'])
    res_blast_sum.index = res_blast_sum.index.set_names(['c_call'])
    res_10x_sum.reset_index(drop = False, inplace = True)
    res_blast_sum.reset_index(drop = False, inplace = True)

    if correct_c_call:
        dat = _correct_c_call(dat, primers_dict=correction_dict)
        res_corrected = pd.DataFrame(dat['c_call'])
        res_corrected_sum = pd.DataFrame(res_corrected['c_call'].value_counts(normalize=True)*100)
        res_corrected_sum['group'] = 'corrected'
        res_corrected_sum.columns = ['counts', 'group']
        res_corrected_sum.index = res_corrected_sum.index.set_names(['c_call'])
        res_corrected_sum.reset_index(drop = False, inplace = True)
        res = pd.concat([res_10x_sum, res_blast_sum, res_corrected_sum])
    else:
        res = pd.concat([res_10x_sum, res_blast_sum])
    res = res.reset_index(drop = True)
    res['c_call'] = res['c_call'].astype('category')
    res['c_call'] = res['c_call'].cat.reorder_categories(sorted(list(set(res['c_call'])), reverse=True))

    dat = _add_cell(dat)
    dat.to_csv(_file, sep = '\t', index=False)
    if plot:
        options.figure_size = figsize
        if correct_c_call:
            p = (ggplot(res, aes(x='c_call', y = 'counts', fill='group'))
                + coord_flip()
                + theme_classic()
                + xlab("c_call")
                + ylab("% c calls")
                + geom_bar(stat="identity", position = 'dodge')
                + scale_fill_manual(values=('#e7e7e7','#86bcb6', '#F28e2b'))
                + theme(legend_title = element_blank()))
        else:
            p = (ggplot(res, aes(x='c_call', y = 'counts', fill='group'))
                + coord_flip()
                + theme_classic()
                + xlab("c_call")
                + ylab("% c calls")
                + geom_bar(stat="identity", position = 'dodge')
                + scale_fill_manual(values=('#e7e7e7','#86bcb6'))
                + theme(legend_title = element_blank()))
        print(p)


def reannotate_genes(data, igblast_db = None, germline = None, org ='human', loci = 'ig', fileformat = 'airr', dirs = None, filtered = False, extended = False, verbose = False, *args):
    """
    Reannotate cellranger fasta files with igblastn and parses to airr/changeo data format.

    Parameters
    ----------
    data : list
        list or sequence of fasta file locations, or folder name containing fasta files. if provided as a single string, it will first be converted to a list; this allows for the function to be run on single/multiple samples.
    igblast_db : str, optional
        path to igblast database folder. Defaults to `$IGDATA` environmental variable.
    germline : str, optional
        path to germline database folder. Defaults to `$GERMLINE` environmental variable.
    org : str
        organism of germline database. Default is 'human'.
    loci : str
        mode for igblastn. Default is 'ig' for BCRs. Also accepts 'tr' for TCRs.
    fileformat: str
        format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.
    dirs : str, optional
        path to input files. will also determine folder structure for outout. Defaults to 'dandelion/data'.
    filtered : bool
        whether or not the to use 'filtered_contig' (True) or 'all_contig' (False) as prefix for output files.
    extended : bool
        whether or not to transfer additional 10X annotions to output file. Default is False.
    verbose :
        whether or not to print the igblast command used in the terminal. Default is False.    
    *args
        passed to `dandelion.preprocessing.ext.assigngenes_igblast` and `dandelion.preprocessing.ext.makedb_igblast`.
    Returns
    ----------
        V(D)J data file in airr/changeo data format.
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
                    raise OSError('Environmental variable GERMLINE must be set. Otherwise, please provide path to folder containing germline fasta files.')
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

def reassign_alleles(data, out_folder, dirs = None, germline = None, org = 'human', fileformat = 'airr', seq_field = 'sequence_alignment', v_field='v_call_genotyped', d_field='d_call', j_field='j_call', germ_types='dmask', plot = True, figsize = (4,3), sample_dict = None, split_write_out = True, filtered = False, out_filename = None, verbose = False):
    """
    Correct allele calls based on a personalized genotype using tigger-reassignAlleles. It uses a subject-specific genotype to correct correct preliminary allele assignments of a set of sequences derived from a single subject.

    Parameters
    ----------
    data : list
        list or sequence of data folders/file locations. if provided as a single string, it will first be converted to a list; this allows for the function to be run on single/multiple samples.
    out_folder : str
        name of folder for concatenated data file and genotyped files.
    dirs : str, optional
        path to input files. will also determine folder structure for outout. Defaults to 'dandelion/data'.
    germline : str, optional
        path to germline database folder. Defaults to `$GERMLINE` environmental variable.
    org : str
        organism of germline database. Default is 'human'.
    fileformat : str
        format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.
    org : str
        organism of germline database. Default is 'human'.
    seq_field : str
        name of column containing the aligned sequence. Default is 'sequence_alignment' (airr).
    v_field : str
        name of column containing the germline V segment call. Default is 'v_call_genotyped' (airr) after tigger.
    d_field : str
        name of column containing the germline d segment call. Default is 'd_call' (airr).
    j_field : str
        name of column containing the germline j segment call. Default is 'j_call' (airr).
    germ_types : str
        Specify type(s) of germlines to include full germline, germline with D segment masked, or germline for V segment only. Default is 'dmask'.
    plot : bool
        whether or not to plot reassignment summary metrics. Default is True.
    figsize : tuple[float, float]
        size of figure. Default is (4, 3).
    sample_dict : dict, optional
        dictionary for creating a sample_id column in the concatenated file.
    split_write_out : bool
        whether or not to write out the processed file. Default is True. If False, it will return a `Dandelion` object instead.
    filtered : bool
        whether or not the to use 'filtered_contig' (True) or 'all_contig' (False) as prefix for output files. Ignored if out_filename is specified.
    out_filename : str, optional
        if provided, will save output to this filename.
    verbose : bool
        Whether or not to print the command used in the terminal. Default is False.
    Returns
    ----------
        Individual V(D)J data files with v_call_genotyped column containing reassigned heavy chain v calls
        Dandelion object holding updated `.data` slot if split_write_out is False.
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

    germline_dict = {'changeo':'_igblast_db-pass_genotype.fasta', 'airr':'_igblast_gap_genotype.fasta'}

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

    # initialise the germline references
    # germline_ref = readGermlines([gml])

    # and now to add it back to the original folders
    sleep(0.5)
    if out_filename is None:
        if filtered:
            out_h = load_data(outDir+'filtered_contig_heavy'+fileformat_dict[fileformat])
            # out = pd.read_csv(outDir+'filtered_contig'+fileformat_dict[fileformat], sep = '\t', dtype = 'object')
            dat_['v_call_genotyped'] = pd.Series(out_h['v_call_genotyped'])
            dat_ = _return_IGKV_IGLV(dat_)
            res = Dandelion(dat_)
            # update with the personalized germline database
            res.update_germline(outDir+'filtered_contig_heavy'+germline_dict[fileformat], germline, org)
            create_germlines(res, germline = germline, org = org, seq_field = seq_field, v_field = v_field, d_field = d_field, j_field = j_field, germ_types = germ_types, fileformat = fileformat)
            print('   Saving corrected genotyped object')
            sleep(0.5)
            res.data.to_csv(outDir+'filtered_contig'+fileformat_dict[fileformat], index = False, sep = '\t')
        else:
            out_h = load_data(outDir+'all_contig_heavy'+fileformat_dict[fileformat])
            # out = pd.read_csv(outDir+'all_contig'+fileformat_dict[fileformat], sep = '\t', dtype = 'object')
            dat_['v_call_genotyped'] = pd.Series(out_h['v_call_genotyped'])
            dat_ = _return_IGKV_IGLV(dat_)
            res = Dandelion(dat_)
            # update with the personalized germline database
            res.update_germline(outDir+'all_contig_heavy'+germline_dict[fileformat], germline, org)
            create_germlines(res, germline = germline, org = org, seq_field = seq_field, v_field = v_field, d_field = d_field, j_field = j_field, germ_types = germ_types, fileformat = fileformat)
            print('   Saving corrected genotyped object')
            sleep(0.5)
            res.data.to_csv(outDir+'all_contig'+fileformat_dict[fileformat], index = False, sep = '\t')
    else:
        out_h = load_data(outDir+'heavy_'+out_filename.replace('.tsv', '_genotyped.tsv'))
        # out = pd.read_csv(outDir+out_filename.replace('.tsv', '_genotyped.tsv'), sep = '\t', dtype = 'object')
        dat_['v_call_genotyped'] = pd.Series(out_h['v_call_genotyped'])
        dat_ = _return_IGKV_IGLV(dat_)
        res = Dandelion(dat_)
        res.update_germline(outDir+'heavy_'+out_filename.replace('.tsv', '.fasta'), germline, org)
        create_germlines(res, germline = germline, org = org, seq_field = seq_field, v_field = v_field, d_field = d_field, j_field = j_field, germ_types = germ_types, fileformat = fileformat)
        print('   Saving corrected genotyped object')
        sleep(0.5)
        res.data.to_csv(out_filename.replace('.tsv', '_genotyped.tsv'), index = False, sep = '\t')

    # reset dat_
    dat_ = res.data.copy()    

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
            res_x = out_h[(out_h['sample_id']==samp)]
            V_ = [re.sub('[*][0-9][0-9]', '', v) for v in res_x['v_call']]
            V_g = [re.sub('[*][0-9][0-9]', '', v) for v in res_x['v_call_genotyped']]
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
        p = (ggplot(final_table, aes(x='sample_id', y = 'var', fill='var_group'))
            + coord_flip()
            + theme_classic()
            + xlab("sample_id")
            + ylab("% allele calls")
            + ggtitle("Genotype reassignment with TIgGER")
            + geom_bar(stat="identity")
            + facet_grid('~'+str('vgroup'), scales="free_y")
            + scale_fill_manual(values=('#86bcb6', '#F28e2b'))
            + theme(legend_title = element_blank()))
        print(p)

    if split_write_out:
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
    else:
        return(res)

def create_germlines(self, germline = None, org = 'human', seq_field='sequence_alignment', v_field='v_call', d_field='d_call', j_field='j_call', germ_types='dmask', fileformat='airr'):
    """
    Runs CreateGermlines.py to reconstruct the germline V(D)J sequence, from which the Ig lineage and mutations can be inferred.

    Parameters
    ----------
    self : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones have been determined.
    germline : str, optional
        path to germline database folder. Defaults to `$GERMLINE` environmental variable.
    org : str
        organism of germline database. Default is 'human'.
    seq_field : str
        name of column containing the aligned sequence. Default is 'sequence_alignment' (airr).
    v_field : str
        name of column containing the germline V segment call. Default is 'v_call' (airr).
    d_field : str
        name of column containing the germline d segment call. Default is 'd_call' (airr).
    j_field : str
        name of column containing the germline j segment call. Default is 'j_call' (airr).
    germ_types : str
        Specify type(s) of germlines to include full germline, germline with D segment masked, or germline for V segment only. Default is 'dmask'.
    fileformat : str
        format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.
    Returns
    ----------
        V(D)J data file with reconstructed germline sequences.
    """
    start = logg.info('Reconstructing germline sequences')
    env = os.environ.copy()
    if germline is None:
        try:
            gml = env['GERMLINE']
        except:
            raise OSError('Environmental variable GERMLINE must be set. Otherwise, please provide path to folder containing germline fasta files.')
        gml = gml+'imgt/'+org+'/vdj/'
    else:
        if os.path.isdir(germline):
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

    def _create_germlines_object(self, references, seq_field, v_field, d_field, j_field, germ_types, fileformat):
        """
        Write germline sequences to tab-delimited database file

        Arguments:
        self : dandelion_class object
        references : folders and/or files containing germline repertoire data in FASTA format.
        seq_field : field in which to look for sequence.
        v_field : field in which to look for V call.
        d_field : field in which to look for D call.
        j_field : field in which to look for J call.
        # cloned : if True build germlines by clone, otherwise build individual germlines.
        # clone_field : field containing clone identifiers; ignored if cloned=False.
        germ_types : list of germline sequence types to be output from the set of 'full', 'dmask', 'vonly', 'regions'
        fileformat : str
            format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.

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

        if type(references) is dict:
            reference_dict = references
        else:
            if type(references) is not list:
                ref = [references]
            else:
                ref = references
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
                # clone_field = schema.toReceptor(clone_field)

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
            # clone_field = schema.toReceptor(clone_field)
            # Define Receptor iterator
            receptor_iter = ((self.loc[x, ].sequence_id, self.loc[x, ]) for x in self.index)
        
        out = {}
        # Iterate over rows
        for key, records in tqdm(receptor_iter, desc = "   Building {} germline sequences".format(germ_types)):
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
            datx = load_data(self.data)
            for x in germline_df.columns:
                datx[x] = pd.Series(germline_df[x])
            self.__init__(data = datx, germline = reference_dict)
        elif self.__class__ == pd.DataFrame:
            datx = load_data(self)
            for x in germline_df.columns:
                datx[x] = pd.Series(germline_df[x])
            output = Dandelion(data = datx, germline = reference_dict)
            return(output)
        logg.info(' finished', time=start,
        deep=('Updated Dandelion object: \n'
        '   \'data\', updated germline alignment in contig-indexed clone table\n'
        '   \'germline\', updated germline reference\n'))


    def _create_germlines_file(file, references, seq_field, v_field, d_field, j_field, germ_types, fileformat):
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
        germ_types : list of germline sequence types to be output from the set of 'full', 'dmask', 'vonly', 'regions'
        fileformat : str
                format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.
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

        if type(references) is dict:
            reference_dict = references
        else:
            if type(references) is not list:
                ref = [references]
            else:
                ref = references
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
        # clone_field = schema.toReceptor(clone_field)
        # Define Receptor iterator
        receptor_iter = ((x.sequence_id, [x]) for x in db_iter)
        
        out = {}
        # Iterate over rows
        for key, records in tqdm(receptor_iter, desc = "   Building {} germline sequences".format(germ_types)):
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

        out = Dandelion(data = file, germline = reference_dict)
        for x in germline_df.columns:
            out.data[x] = pd.Series(germline_df[x])

        if os.path.isfile(str(file)):
            out.data.to_csv("{}/{}_germline_{}.tsv".format(os.path.dirname(file), os.path.basename(file).split('.tsv')[0], germ_types), sep = '\t', index = False)        
        return(out)

    if (type(germline) is dict) or (type(germline) is list):
        if self.__class__ == Dandelion:
            _create_germlines_object(self, germline, seq_field, v_field, d_field, j_field, germ_types, fileformat)
        elif self.__class__ == pd.DataFrame:
            return(_create_germlines_object(self, germline, seq_field, v_field, d_field, j_field, germ_types, fileformat))
        else:
            return(_create_germlines_file(self, germline, seq_field, v_field, d_field, j_field, germ_types, fileformat))
    else:
        if self.__class__ == Dandelion:
            if len(self.germline) is not 0:
                _create_germlines_object(self, self.germline, seq_field, v_field, d_field, j_field, germ_types, fileformat)
            else:
                _create_germlines_object(self, gml, seq_field, v_field, d_field, j_field, germ_types, fileformat)
        elif self.__class__ == pd.DataFrame:
            return(_create_germlines_object(self, gml, seq_field, v_field, d_field, j_field, germ_types, fileformat))
        else:
            return(_create_germlines_file(self, gml, seq_field, v_field, d_field, j_field, germ_types, fileformat))

def recipe_scanpy_qc(self, max_genes=2500, min_genes=200, mito_cutoff=0.05, pval_cutoff=0.1, min_counts=None, max_counts=None):
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
        Maximum percentage mitochondrial content allowed for a cell to pass filtering. Default is 0.05.
    pval_cutoff : float
        Maximum Benjamini-Hochberg corrected p value from doublet detection protocol allowed for a cell to pass filtering. Default is 0.05.
    min_counts : int
        Minimum number of counts required for a cell to pass filtering. Default is None.
    max_counts : int
        Maximum number of counts required for a cell to pass filtering. Default is None.
    Returns
    -------
        `AnnData` of shape n_obs × n_vars where obs now contain filtering information. Rows correspond to cells and columns to genes.

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

def filter_bcr(data, adata, filter_bcr=True, filter_rna=True, rescue_igh=True, umi_foldchange_cutoff=2, filter_lightchains=True, filter_missing=True, outdir=None, outFilePrefix=None, filtered=False):
    """
    Filters doublets and poor quality cells and corresponding contigs based on provided V(D)J `DataFrame` and `AnnData` objects. Depends on a `AnnData`.obs slot populated with 'filter_rna' column.
    If the aligned sequence is an exact match between contigs, the contigs will be merged into the one with the highest umi count, adding the summing the umi count of the duplicated contigs to duplicate_count column. After this check, if there are still multiple contigs, cells with multiple IGH contigs are filtered unless `rescue_igh` is True, where by the umi counts for each IGH contig will then be compared. The contig with the highest umi that is > umi_foldchange_cutoff (default is empirically set at 5) from the lowest will be retained.
    If there's multiple contigs that survive the 'rescue', then all contigs will be filtered. The default behaviour is to also filter cells with multiple lightchains but this may sometimes be a true biological occurrence; toggling filter_lightchains to False will rescue the mutltiplet light chains.
    Lastly, contigs with no corresponding cell barcode in the AnnData object is filtered if filter_missing is True. However, this may be useful to toggle to False if more contigs are preferred to be kept or for integrating with bulk reperotire seq data.

    Parameters
    ----------
    data : DataDrame, str
        V(D)J airr/changeo data to filter. Can be pandas `DataFrame` object or file path as string.
    adata : AnnData
        AnnData object to filter.
    filter_bcr : bool
        If True, V(D)J `DataFrame` object returned will be filtered. Default is True.
    filter_rna : bool
        If True, `AnnData` object returned will be filtered. Default is True.
    rescue_igh : bool
        If True, rescues IGH contigs with highest umi counts with a requirement that it passes the `umi_foldchange_cutoff` option. In addition, the sum of the all the heavy chain contigs must be greater than 3 umi or all contigs will be filtered. Default is True.
    umi_foldchange_cutoff : int
        related to minimum fold change required to rescue heavy chain contigs/barcode otherwise they will be marked as doublets. Default is empirically set at 2-fold. 
    filter_lightchains : bool
        cells with multiple light chains will be marked to filter. Default is True.
    filter_missing : bool
        cells in V(D)J data not found in `AnnData` object will be marked to filter. Default is True. This may be useful for toggling to False if integrating with bulk data.
    outdir : str, optional
        If specified, out file will be in this location
    outFilePrefix : str, optional
        If specified, the out file name will have this prefix
    filtered : bool
        If True, will create filenames with 'filtered_contig' as prefix. if False, will create filenames with 'all_contig' as prefix. ignored if outFilePrefix is specified.
    Returns
    -------
        V(D)J `DataFrame` object in airr/changeo format and `AnnData` object.
    """
    dat = load_data(data)
    h = Tree()
    l = Tree()
    h_umi = Tree()
    h_dup = Tree()
    l_umi = Tree()
    h_seq = Tree()
    l_seq = Tree()
    poor_qual, h_doublet, l_doublet, drop_contig  = [], [], [], []
    
    locus_dict = dict(zip(dat['sequence_id'],dat['locus']))
    barcode = list(set(dat['cell_id']))

    if 'filter_rna' not in adata.obs:
        raise TypeError("AnnData obs does not contain 'filter_rna' column. Please run `pp.recipe_scanpy_qc` before continuing.")

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

    # rather than leaving a nan cell, i will create a 0 column for now
    dat['duplicate_count'] = 0
    for b in tqdm(barcode, desc = 'Marking barcodes with poor quality BCRs and BCR doublets'):
        hc_id = list(dat[(dat['cell_id'].isin([b])) & (dat['locus'] == 'IGH')]['sequence_id'])        
        hc_umi = [int(x) for x in dat[(dat['cell_id'].isin([b])) & (dat['locus'] == 'IGH')]['umi_count']]
        hc_seq = [x for x in dat[(dat['cell_id'].isin([b])) & (dat['locus'] == 'IGH')]['sequence_alignment']]
        hc_dup = [int(x) for x in dat[(dat['cell_id'].isin([b])) & (dat['locus'] == 'IGH')]['duplicate_count']]

        lc_id = list(dat[(dat['cell_id'].isin([b])) & (dat['locus'].isin(['IGK', 'IGL']))]['sequence_id'])
        lc_umi = [int(x) for x in dat[(dat['cell_id'].isin([b])) & (dat['locus'].isin(['IGK', 'IGL']))]['umi_count']]
        lc_seq = [x for x in dat[(dat['cell_id'].isin([b])) & (dat['locus'].isin(['IGK', 'IGL']))]['sequence_alignment']]
        
        h[b] = hc_id
        h_umi[b] = hc_umi
        h_seq[b] = hc_seq
        h_dup[b] = hc_dup

        l[b] = lc_id
        l_umi[b] = lc_umi
        l_seq[b] = lc_seq
        
        # marking doublets defined by heavy chains
        if len(h[b]) > 1:
            if len(list(set(h_seq[b]))) == 1:          
                highest_umi_h = max(h_umi[b])
                highest_umi_h_idx = [i for i, j in enumerate(h_umi[b]) if j == highest_umi_h]
                drop_contig.append(h[b][~highest_umi_h_idx[0]])                    
                keep_hc_contig = h[b][highest_umi_h_idx[0]]
                dat.loc[keep_hc_contig, 'duplicate_count'] = int(np.sum(h_umi[b][~highest_umi_h_idx[0]]))
                    
                hc_id = list(dat[(dat['cell_id'].isin([b])) & (dat['locus'] == 'IGH')]['sequence_id'])
                hc_umi = [int(x) for x in dat[(dat['cell_id'].isin([b])) & (dat['locus'] == 'IGH')]['umi_count']]
                hc_dup = [int(x) for x in dat[(dat['cell_id'].isin([b])) & (dat['locus'] == 'IGH')]['duplicate_count']]
                h[b] = hc_id
                h_umi[b] = hc_umi
                h_dup[b] = hc_dup
                h_seq[b] = hc_seq
            if len(h[b]) > 1:
                if rescue_igh:
                    highest_umi = max(h_umi[b])
                    lowest_umi = min(h_umi[b])
                    sum_umi = sum(h_umi[b]+h_dup[b])
                    highest_umi_idx = [i for i, j in enumerate(h_umi[b]) if j == highest_umi]
                    if len(highest_umi_idx) > 1:
                        h_doublet.append(b)
                    if sum_umi < 4:
                        h_doublet.append(b)
                    if highest_umi/lowest_umi < umi_foldchange_cutoff:
                        h_doublet.append(b)
                    if len(highest_umi_idx) == 1 and highest_umi/lowest_umi >= umi_foldchange_cutoff:
                        drop_contig.append(h[b][~highest_umi_idx[0]])
                else:
                    h_doublet.append(b)

        if len(l[b]) > 1:
            if len(list(set(l_seq[b]))) == 1:
                highest_umi_l = max(l_umi[b])
                highest_umi_l_idx = [i for i, j in enumerate(l_umi[b]) if j == highest_umi_l]
                drop_contig.append(l[b][~highest_umi_l_idx[0]])
                keep_lc_contig = l[b][highest_umi_l_idx[0]]
                dat.loc[keep_lc_contig, 'duplicate_count'] = int(np.sum(l_umi[b][~highest_umi_l_idx[0]]))
                lc_id = list(dat[(dat['cell_id'].isin([b])) & (dat['locus'].isin(['IGK', 'IGL']))]['sequence_id'])
                lc_umi = [int(x) for x in dat[(dat['cell_id'].isin([b])) & (dat['locus'].isin(['IGK', 'IGL']))]['umi_count']]                
                l[b] = lc_id
                l_umi[b] = lc_umi
                l_seq[b] = lc_seq
        
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
        _dat = _dat[~(_dat['sequence_id'].isin(drop_contig))]

        if os.path.isfile(str(data)):
            _dat.to_csv("{}/{}_filtered.tsv".format(os.path.dirname(data), os.path.basename(data).split('.tsv')[0]), sep = '\t', index = None)
        else:
            if filtered:
                outFile_prefix = 'filtered_contig'
            else:
                outFile_prefix = 'all_contig'

            if (outdir is None) & (outFilePrefix is not None):
                _dat.to_csv("{}/{}_filtered.tsv".format('dandelion/data', str(outFilePrefix)), sep = '\t', index = None)
            elif (outdir is not None) & (outFilePrefix is None):
                _dat.to_csv("{}/{}_filtered.tsv".format(str(outdir), outFile_prefix), sep = '\t', index = None)
            elif (outdir is None) & (outFilePrefix is None):
                _dat.to_csv("{}/{}_filtered.tsv".format('dandelion/data', outFile_prefix), sep = '\t', index = None)
            elif (outdir is not None) & (outFilePrefix is None):
                _dat.to_csv("{}/{}_filtered.tsv".format(str(outdir), outFile_prefix), sep = '\t', index = None)

    if filter_rna:
        _adata = adata[~(adata.obs_names.isin(filter_ids))] # not saving the scanpy object because there's no need to at the moment
    else:
        _adata = adata.copy()

    return(_dat, _adata)