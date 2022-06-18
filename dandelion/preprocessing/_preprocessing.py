#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 17:56:02
# @Last Modified by:   Kelvin
# @Last Modified time: 2022-06-18 14:31:22
"""preprocessing module."""
import anndata as ad
import functools
import numpy as np
import os
import pandas as pd
import re

from Bio import Align
from changeo.Gene import buildGermline
from changeo.IO import getFormatOperators, readGermlines, checkFields
from changeo.Receptor import AIRRSchema, ChangeoSchema, Receptor, ReceptorData
from collections import OrderedDict
from os import PathLike
from plotnine import (
    ggplot,
    geom_bar,
    geom_col,
    ggtitle,
    scale_fill_manual,
    coord_flip,
    options,
    element_blank,
    aes,
    xlab,
    ylab,
    facet_wrap,
    facet_grid,
    theme_classic,
    theme,
    annotate,
    theme_bw,
    geom_histogram,
    geom_vline,
    save_as_pdf_pages,
)
from scanpy import logging as logg
from subprocess import run
from time import sleep
from tqdm import tqdm
from typing import Union, Sequence, Tuple, Optional

from .external._preprocessing import (
    assigngenes_igblast,
    makedb_igblast,
    parsedb_heavy,
    parsedb_light,
    tigger_genotype,
    creategermlines,
)
from ..utilities._core import *
from ..utilities._io import *
from ..utilities._utilities import *


TRUES = ["T", "True", "true", "TRUE", True]
FALSES = ["F", "False", "false", "FALSE", False]
HEAVYLONG = ["IGH", "TRB", "TRD"]
LIGHTSHORT = ["IGK", "IGL", "TRA", "TRG"]


def format_fasta(
    fasta: Union[PathLike, str],
    prefix: Optional[str] = None,
    suffix: Optional[str] = None,
    sep: Optional[str] = None,
    remove_trailing_hyphen_number: bool = True,
    high_confidence_filtering: bool = False,
    outdir: Optional[str] = None,
    filename_prefix: Optional[str] = None,
):
    """
    Add prefix to the headers/contig ids in input fasta and annotation file.

    Parameters
    ----------
    fasta : str
        path to fasta file.
    prefix : str, Optional
        prefix to append to the headers/contig ids.
    suffix : str, Optional
        suffix to append to the headers/contig ids.
    sep : str, Optional
        separator after prefix or before suffix to append to the headers/contig
        ids.
    remove_trailing_hyphen_number : bool
        whether or not to remove the trailing hyphen number e.g. '-1' from the
        ell/contig barcodes.
    high_confidence_filtering : bool
        whether ot not to filter to only `high confidence` contigs.
    outdir : str, Optional
        path to output location. None defaults to 'dandelion'.
    filename_prefix : str, Optional
        prefix of file name preceding '_contig'. None defaults to 'filtered'.

    Returns
    -------
    Formatted fasta file with new headers containing prefix
    """
    if filename_prefix is None:
        filename_pre = "filtered"
    else:
        filename_pre = filename_prefix

    filePath = None
    filePath = check_fastapath(fasta, filename_prefix=filename_pre)

    if filePath is None:
        raise FileNotFoundError(
            "Path to fasta file is unknown. Please "
            + "specify path to fasta file or folder containing fasta file. "
            + "Starting folder should only contain 1 fasta file."
        )

    fh = open(filePath, "r")
    seqs = {}

    if sep is None:
        separator = "_"
    else:
        separator = str(sep)

    for header, sequence in fasta_iterator(fh):
        if prefix is None and suffix is None:
            seqs[header] = sequence
        elif prefix is not None:
            if suffix is not None:
                if remove_trailing_hyphen_number:
                    newheader = (
                        str(prefix)
                        + separator
                        + str(header).split("_contig")[0].split("-")[0]
                        + separator
                        + str(suffix)
                        + "_contig"
                        + str(header).split("_contig")[1]
                    )
                else:
                    newheader = (
                        str(prefix)
                        + separator
                        + str(header).split("_contig")[0]
                        + separator
                        + str(suffix)
                        + "_contig"
                        + str(header).split("_contig")[1]
                    )
            else:
                if remove_trailing_hyphen_number:
                    newheader = (
                        str(prefix)
                        + separator
                        + str(header).split("_contig")[0].split("-")[0]
                        + "_contig"
                        + str(header).split("_contig")[1]
                    )
                else:
                    newheader = str(prefix) + separator + str(header)
            seqs[newheader] = sequence
        else:
            if suffix is not None:
                if remove_trailing_hyphen_number:
                    newheader = (
                        str(header).split("_contig")[0].split("-")[0]
                        + separator
                        + str(suffix)
                        + "_contig"
                        + str(header).split("_contig")[1]
                    )
                else:
                    newheader = (
                        str(header).split("_contig")[0]
                        + separator
                        + str(suffix)
                        + "_contig"
                        + str(header).split("_contig")[1]
                    )
            else:
                newheader = str(header)
            seqs[newheader] = sequence

    fh.close()

    if os.path.isfile(filePath):
        basedir = os.path.dirname(filePath)
    elif os.path.isdir(filePath):
        basedir = os.path.dirname(filePath)
    else:
        basedir = os.getcwd()

    if outdir is None:
        out_dir = basedir.rstrip("/") + "/" + "dandelion/"
    else:
        if not outdir.endswith("/"):
            out_dir = outdir + "/"

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # format the barcode and contig_id in the corresponding annotation file too
    anno = (
        basedir
        + "/"
        + os.path.basename(filePath).replace(".fasta", "_annotations.csv")
    )
    data = pd.read_csv(anno, dtype="object")

    if prefix is not None:
        if suffix is not None:
            if remove_trailing_hyphen_number:
                data["contig_id"] = [
                    str(prefix)
                    + separator
                    + str(c).split("_contig")[0].split("-")[0]
                    + separator
                    + str(suffix)
                    + "_contig"
                    + str(c).split("_contig")[1]
                    for c in data["contig_id"]
                ]
                data["barcode"] = [
                    str(prefix)
                    + separator
                    + str(b).split("-")[0]
                    + separator
                    + str(suffix)
                    for b in data["barcode"]
                ]
            else:
                data["contig_id"] = [
                    str(prefix)
                    + separator
                    + str(c).split("_contig")[0]
                    + separator
                    + str(suffix)
                    + "_contig"
                    + str(c).split("_contig")[1]
                    for c in data["contig_id"]
                ]
                data["barcode"] = [
                    str(prefix) + separator + str(b) + separator + str(suffix)
                    for b in data["barcode"]
                ]
        else:
            if remove_trailing_hyphen_number:
                data["contig_id"] = [
                    str(prefix)
                    + separator
                    + str(c).split("_contig")[0].split("-")[0]
                    + "_contig"
                    + str(c).split("_contig")[1]
                    for c in data["contig_id"]
                ]
                data["barcode"] = [
                    str(prefix) + separator + str(b).split("-")[0]
                    for b in data["barcode"]
                ]
            else:
                data["contig_id"] = [
                    str(prefix) + separator + str(c) for c in data["contig_id"]
                ]
                data["barcode"] = [
                    str(prefix) + separator + str(b) for b in data["barcode"]
                ]
    else:
        if suffix is not None:
            if remove_trailing_hyphen_number:
                data["contig_id"] = [
                    str(c).split("_contig")[0].split("-")[0]
                    + separator
                    + str(suffix)
                    + "_contig"
                    + str(c).split("_contig")[1]
                    for c in data["contig_id"]
                ]
                data["barcode"] = [
                    str(b).split("-")[0] + separator + str(suffix)
                    for b in data["barcode"]
                ]
            else:
                data["contig_id"] = [
                    str(c).split("_contig")[0]
                    + separator
                    + str(suffix)
                    + "_contig"
                    + str(c).split("_contig")[1]
                    for c in data["contig_id"]
                ]
                data["barcode"] = [
                    str(b) + separator + str(suffix) for b in data["barcode"]
                ]
        else:
            data["contig_id"] = [str(c) for c in data["contig_id"]]
            data["barcode"] = [str(b) for b in data["barcode"]]

    out_anno = out_dir + os.path.basename(filePath).replace(
        ".fasta", "_annotations.csv"
    )
    out_fasta = out_dir + os.path.basename(filePath)
    fh1 = open(out_fasta, "w")
    fh1.close()
    out = ""

    if high_confidence_filtering:
        hiconf_contigs = [
            x
            for x, y in zip(data["contig_id"], data["high_confidence"])
            if y in TRUES
        ]
        seqs = {hiconf: seqs[hiconf] for hiconf in hiconf_contigs}

        data = data[data["contig_id"].isin(hiconf_contigs)]

    for l in seqs:
        out = ">" + l + "\n" + seqs[l] + "\n"
        Write_output(out, out_fasta)
    data.to_csv(out_anno, index=False)


def format_fastas(
    fastas: Sequence,
    prefix: Optional[Sequence] = None,
    suffix: Optional[Sequence] = None,
    sep: Optional[str] = None,
    remove_trailing_hyphen_number: bool = True,
    high_confidence_filtering: bool = False,
    outdir: Optional[str] = None,
    filename_prefix: Optional[Union[Sequence, str]] = None,
):
    """
    Add prefix to the headers/contig ids in input fasta and annotation file.

    Parameters
    ----------
    fastas : Sequence
        list of paths to fasta files.
    prefix : list, Optional
        list of prefixes to append to headers/contig ids in each fasta file.
    suffix : str, Optional
        list of suffixes to append to headers/contig ids in each fasta file.
    sep : str, Optional
        separator after prefix or before suffix to append to the headers/contig
        ids.
    remove_trailing_hyphen_number : bool
        whether or not to remove the trailing hyphen number e.g. '-1' from the
        cell/contig barcodes.
    high_confidence_filtering : bool
        whether ot not to filter to only `high confidence` contigs.
    outdir : str, Optional
        path to out put location. Default is None, which is 'dandelion'.
    filename_prefix : str, Optional
        list of prefixes of file names preceding '_contig'. None defaults to
        'filtered'.

    Returns
    -------
    Formatted fasta file with new headers containing prefix
    """
    if type(fastas) is not list:
        fastas = [fastas]
    if type(filename_prefix) is not list:
        filename_prefix = [filename_prefix]
    if all(t is None for t in filename_prefix):
        filename_prefix = [None for f in fastas]

    if prefix is not None:
        if type(prefix) is not list:
            prefix = [prefix]
        prefix_dict = dict(zip(fastas, prefix))
    if suffix is not None:
        if type(suffix) is not list:
            suffix = [suffix]
        suffix_dict = dict(zip(fastas, suffix))

    for i in tqdm(range(0, len(fastas)), desc="Formating fasta(s) "):
        if prefix is None and suffix is None:
            format_fasta(
                fastas[i],
                prefix=None,
                suffix=None,
                sep=None,
                remove_trailing_hyphen_number=remove_trailing_hyphen_number,
                high_confidence_filtering=high_confidence_filtering,
                outdir=outdir,
                filename_prefix=filename_prefix[i],
            )
        elif prefix is not None:
            if suffix is not None:
                format_fasta(
                    fastas[i],
                    prefix=prefix_dict[fastas[i]],
                    suffix=suffix_dict[fastas[i]],
                    sep=sep,
                    remove_trailing_hyphen_number=remove_trailing_hyphen_number,
                    high_confidence_filtering=high_confidence_filtering,
                    outdir=outdir,
                    filename_prefix=filename_prefix[i],
                )
            else:
                format_fasta(
                    fastas[i],
                    prefix=prefix_dict[fastas[i]],
                    suffix=None,
                    sep=sep,
                    remove_trailing_hyphen_number=remove_trailing_hyphen_number,
                    high_confidence_filtering=high_confidence_filtering,
                    outdir=outdir,
                    filename_prefix=filename_prefix[i],
                )
        else:
            if suffix is not None:
                format_fasta(
                    fastas[i],
                    prefix=None,
                    suffix=suffix_dict[fastas[i]],
                    sep=sep,
                    remove_trailing_hyphen_number=remove_trailing_hyphen_number,
                    high_confidence_filtering=high_confidence_filtering,
                    outdir=outdir,
                    filename_prefix=filename_prefix[i],
                )
            else:
                format_fasta(
                    fastas[i],
                    prefix=None,
                    suffix=None,
                    sep=None,
                    remove_trailing_hyphen_number=remove_trailing_hyphen_number,
                    high_confidence_filtering=high_confidence_filtering,
                    outdir=outdir,
                    filename_prefix=filename_prefix[i],
                )


def assign_isotype(
    fasta: Union[str, PathLike],
    fileformat: Literal["blast", "changeo", "airr"] = "blast",
    org: Literal["human", "mouse"] = "human",
    evalue: float = 1e-4,
    correct_c_call: bool = True,
    correction_dict: Union[Dict, None] = None,
    plot: bool = True,
    save_plot: bool = False,
    show_plot: bool = True,
    figsize: Tuple[Union[int, float], Union[int, float]] = (4, 4),
    blastdb: Optional[str] = None,
    allele: bool = False,
    filename_prefix: Optional[str] = None,
    verbose: bool = False,
):
    """
    Annotate contigs with constant region call using blastn.

    Parameters
    ----------
    fasta : str, PathLike
        path to fasta file.
    fileformat : str
        format of V(D)J file/objects. Default is 'blast'. Also accepts
        'changeo' (same behaviour as 'blast') and 'airr'.
    org : str
        organism of reference folder. Default is 'human'.
    evalue : float
        This is the statistical significance threshold for reporting matches
        against database sequences. Lower EXPECT thresholds are more stringent
        and report only high similarity matches. Choose higher EXPECT value
        (for example 1 or more) if you expect a low identity between your query
        sequence and the targets.
    correct_c_call : bool
        whether or not to adjust the c_calls after blast based on provided
        primers specified in `primer_dict` option. Default is True.
    correction_dict : Dict, Optional
        a nested dictionary contain isotype/c_genes as keys and primer
        sequences as records to use for correcting annotated c_calls. Defaults
        to a curated dictionary for human sequences if left as none.
    plot : bool
        whether or not to plot reassignment summary metrics. Default is True.
    save_plot : bool
        whether or not to save plot.
    show_plot : bool
        whether or not to show plot.
    figsize : Tuple[Union[int,float], Union[int,float]]
        size of figure. Default is (4, 4).
    blastdb : str, Optional
        path to blast database. Defaults to `$BLASTDB` environmental variable.
    allele : bool
        whether or not to return allele calls. Default is False.
    filename_prefix : str, Optional
        prefix of file name preceding '_contig'. None defaults to 'filtered'.
    verbose : bool
        whether or not to print the blast command in terminal.
        Default is False.

    Returns
    -------
    V(D)J tsv files with constant genes annotated.
    """
    aligner = Align.PairwiseAligner()

    def two_gene_correction(self, i, dictionary):
        """Pairwise alignmet for two genes."""
        key1, key2 = dictionary.keys()
        seq = self.loc[i, "c_sequence_alignment"].replace("-", "")
        alignments1 = aligner.align(dictionary[key1], seq)
        alignments2 = aligner.align(dictionary[key2], seq)
        score1 = alignments1.score
        score2 = alignments2.score
        if score1 == score2:
            self.at[i, "c_call"] = str(key1) + "," + str(key2)
        if score1 > score2:
            self.at[i, "c_call"] = str(key1)
        if score1 < score2:
            self.at[i, "c_call"] = str(key2)

    def three_gene_correction(self, i, dictionary):
        """Pairwise alignmet for three genes."""
        key1, key2, key3 = dictionary.keys()
        seq = self.loc[i, "c_sequence_alignment"].replace("-", "")
        alignments1 = aligner.align(dictionary[key1], seq)
        alignments2 = aligner.align(dictionary[key2], seq)
        alignments3 = aligner.align(dictionary[key3], seq)
        score1 = alignments1.score
        score2 = alignments2.score
        score3 = alignments3.score
        if score1 == score2 == score3:
            self.at[i, "c_call"] = str(key1) + "," + str(key2) + "," + str(key3)
        elif score1 > score2 and score1 > score3:
            self.at[i, "c_call"] = str(key1)
        elif score2 > score1 and score2 > score3:
            self.at[i, "c_call"] = str(key2)
        elif score3 > score1 and score3 > score2:
            self.at[i, "c_call"] = str(key3)
        elif score1 == score2 and score1 > score3:
            self.at[i, "c_call"] = str(key1) + "," + str(key2)
        elif score1 > score2 and score1 == score3:
            self.at[i, "c_call"] = str(key1) + "," + str(key3)
        elif score2 > score1 and score2 == score3:
            self.at[i, "c_call"] = str(key2) + "," + str(key3)

    def four_gene_correction(self, i, dictionary):
        """Pairwise alignmet for four genes."""
        key1, key2, key3, key4 = dictionary.keys()
        seq = self.loc[i, "c_sequence_alignment"].replace("-", "")
        alignments1 = aligner.align(dictionary[key1], seq)
        alignments2 = aligner.align(dictionary[key2], seq)
        alignments3 = aligner.align(dictionary[key3], seq)
        alignments4 = aligner.align(dictionary[key4], seq)
        score1 = alignments1.score
        score2 = alignments2.score
        score3 = alignments3.score
        score4 = alignments4.score
        if score1 == score2 == score3 == score4:
            self.at[i, "c_call"] = (
                str(key1) + "," + str(key2) + "," + str(key3) + "," + str(key4)
            )
        elif score1 > score2 and score1 > score3 and score1 > score4:
            self.at[i, "c_call"] = str(key1)
        elif score2 > score1 and score2 > score3 and score2 > score4:
            self.at[i, "c_call"] = str(key2)
        elif score3 > score1 and score3 > score2 and score3 > score4:
            self.at[i, "c_call"] = str(key3)
        elif score4 > score1 and score4 > score2 and score4 > score3:
            self.at[i, "c_call"] = str(key4)
        elif score1 == score2 and score1 > score3 and score1 > score4:
            self.at[i, "c_call"] = str(key1) + "," + str(key2)
        elif score1 > score2 and score1 == score3 and score1 > score4:
            self.at[i, "c_call"] = str(key1) + "," + str(key3)
        elif score1 > score2 and score1 > score3 and score1 == score4:
            self.at[i, "c_call"] = str(key1) + "," + str(key4)
        elif score2 == score3 and score2 > score1 and score2 > score4:
            self.at[i, "c_call"] = str(key1) + "," + str(key3)
        elif score2 == score4 and score2 > score1 and score2 > score3:
            self.at[i, "c_call"] = str(key2) + "," + str(key4)
        elif score3 == score4 and score3 > score1 and score3 > score2:
            self.at[i, "c_call"] = str(key3) + "," + str(key4)
        elif score1 == score2 == score3 and score1 > score4:
            self.at[i, "c_call"] = str(key1) + "," + str(key2) + "," + str(key3)
        elif score1 == score2 == score4 and score1 > score3:
            self.at[i, "c_call"] = str(key1) + "," + str(key2) + "," + str(key4)
        elif score1 == score3 == score4 and score1 > score2:
            self.at[i, "c_call"] = str(key1) + "," + str(key3) + "," + str(key4)
        elif score2 == score3 == score4 and score2 > score1:
            self.at[i, "c_call"] = str(key2) + "," + str(key3) + "," + str(key4)

    def _correct_c_call(data, primers_dict=None):
        """Pairiwise alignment for c genes."""
        dat = data.copy()
        if primers_dict is None:
            primer_dict = {
                "IGHG": {
                    "IGHG1": "GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCCTGGGC",
                    "IGHG2": "GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGC",
                    "IGHG3": "GCTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGGGCACAGCGGCCCTGGGC",
                    "IGHG4": "GCTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGGC",
                },
                "IGHA": {
                    "IGHA1": "GCATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCAGCACCCAGCCAGATGGGAACGTGGTCATCGCCTGC",
                    "IGHA2": "GCATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACAGCACCCCCCAAGATGGGAACGTGGTCGTCGCATGC",
                },
                "IGLC7": {
                    "IGLC": "GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAA",
                    "IGLC7": "GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCACCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCGTAA",
                },
                "IGLC3": {
                    "IGLC": "GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAA",
                    "IGLC3": "GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCACCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAA",
                },
                "IGLC6": {
                    "IGLC": "TCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCA",
                    "IGLC6": "TCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGCCTGA",
                },
            }
        else:
            primer_dict = primers_dict

        for i in dat.index:
            if (dat.loc[i, "c_call"] == dat.loc[i, "c_call"]) & (
                dat.loc[i, "c_call"] is not None
            ):
                for k in primer_dict:
                    if k in dat.loc[i, "c_call"]:
                        if len(primer_dict[k]) == 2:
                            two_gene_correction(dat, i, primer_dict[k])
                        elif len(primer_dict[k]) == 3:
                            three_gene_correction(dat, i, primer_dict[k])
                        elif len(primer_dict[k]) == 4:
                            four_gene_correction(dat, i, primer_dict[k])
        return dat

    # main function from here
    format_dict = {
        "changeo": "_igblast_db-pass",
        "blast": "_igblast_db-pass",
        "airr": "_igblast_gap",
    }

    filePath = check_filepath(
        fasta, filename_prefix=filename_prefix, endswith=".fasta"
    )
    if filePath is None:
        raise FileNotFoundError(
            (
                "Path to fasta file is unknown. Please specify path to "
                + "fasta file or folder containing fasta file."
            )
        )

    blast_out = run_blastn(
        fasta=filePath,
        database=blastdb,
        org=org,
        loci="ig",
        call="c",
        max_hsps=1,
        evalue=evalue,
        outfmt=(
            "6 qseqid sseqid pident length mismatch gapopen "
            + "qstart qend sstart send evalue bitscore qseq sseq"
        ),
        dust="no",
        verbose=verbose,
    )
    blast_out.drop_duplicates(subset="sequence_id", keep="first", inplace=True)

    _file = "{}/tmp/{}_genotyped.tsv".format(
        os.path.dirname(filePath),
        os.path.basename(filePath).split(".fasta")[0] + format_dict[fileformat],
    )
    _airrfile = "{}/tmp/{}.tsv".format(
        os.path.dirname(filePath),
        os.path.basename(filePath).split(".fasta")[0] + "_igblast",
    )
    _file2 = "{}/{}_genotyped.tsv".format(
        os.path.dirname(filePath),
        os.path.basename(filePath).split(".fasta")[0] + format_dict[fileformat],
    )

    if verbose:
        print("Loading 10X annotations \n")
    try:
        dat_10x = load_data(_file)
    except FileNotFoundError:
        # maybe a short cut to skip reassign_alleles?
        _file = "{}/tmp/{}.tsv".format(
            os.path.dirname(filePath),
            os.path.basename(filePath).split(".fasta")[0]
            + format_dict[fileformat],
        )
        dat_10x = load_data(_file)
    res_10x = pd.DataFrame(dat_10x["c_call"])
    res_10x["c_call"] = res_10x["c_call"].fillna(value="None")
    if verbose:
        print("Preparing new calls \n")
    dat = load_data(_file)
    for col in [
        "c_call",
        "c_sequence_alignment",
        "c_germline_alignment",
        "c_sequence_start",
        "c_sequence_end",
        "c_score",
        "c_identity",
    ]:
        dat[col] = pd.Series(blast_out[col])
    res_blast = pd.DataFrame(dat["c_call"])
    res_blast = res_blast.fillna(value="None")

    res_10x_sum = pd.DataFrame(
        res_10x["c_call"].value_counts(normalize=True) * 100
    )
    res_blast_sum = pd.DataFrame(
        res_blast["c_call"].value_counts(normalize=True) * 100
    )
    res_10x_sum["group"] = "10X"
    res_blast_sum["group"] = "blast"
    res_10x_sum.columns = ["counts", "group"]
    res_blast_sum.columns = ["counts", "group"]
    res_10x_sum.index = res_10x_sum.index.set_names(["c_call"])
    res_blast_sum.index = res_blast_sum.index.set_names(["c_call"])
    res_10x_sum.reset_index(drop=False, inplace=True)
    res_blast_sum.reset_index(drop=False, inplace=True)

    if (
        correct_c_call
    ):  # TODO: figure out if i need to set up a None correction?
        if verbose:
            print("Correcting C calls \n")
        dat = _correct_c_call(dat, primers_dict=correction_dict)
        res_corrected = pd.DataFrame(dat["c_call"])
        res_corrected = res_corrected.fillna(value="None")
        res_corrected_sum = pd.DataFrame(
            res_corrected["c_call"].value_counts(normalize=True) * 100
        )
        res_corrected_sum["group"] = "corrected"
        res_corrected_sum.columns = ["counts", "group"]
        res_corrected_sum.index = res_corrected_sum.index.set_names(["c_call"])
        res_corrected_sum.reset_index(drop=False, inplace=True)
        res = pd.concat([res_10x_sum, res_blast_sum, res_corrected_sum])
    else:
        res = pd.concat([res_10x_sum, res_blast_sum])

    res = res.reset_index(drop=True)
    res["c_call"] = res["c_call"].fillna(value="None")
    res["c_call"] = [re.sub("[*][0-9][0-9]", "", c) for c in res["c_call"]]
    res["c_call"] = res["c_call"].astype("category")
    res["c_call"] = res["c_call"].cat.reorder_categories(
        sorted(list(set(res["c_call"])), reverse=True)
    )

    if verbose:
        print("Finishing up \n")
    dat["c_call_10x"] = pd.Series(dat_10x["c_call"])
    # some minor adjustment to the final output table
    airr_output = load_data(_airrfile)
    cols_to_merge = [
        "junction_aa_length",
        "fwr1_aa",
        "fwr2_aa",
        "fwr3_aa",
        "fwr4_aa",
        "cdr1_aa",
        "cdr2_aa",
        "cdr3_aa",
        "sequence_alignment_aa",
        "v_sequence_alignment_aa",
        "d_sequence_alignment_aa",
        "j_sequence_alignment_aa",
    ]
    for x in cols_to_merge:
        dat[x] = pd.Series(airr_output[x])

    # remove allellic calls
    dat["c_call"] = dat["c_call"].fillna(value="")
    dat["c_call"] = [re.sub("[*][0-9][0-9]", "", c) for c in dat["c_call"]]

    write_airr(dat, _file2)
    if plot:
        options.figure_size = figsize
        if correct_c_call:
            p = (
                ggplot(res, aes(x="c_call", y="counts", fill="group"))
                + coord_flip()
                + theme_classic()
                + xlab("c_call")
                + ylab("% c calls")
                + geom_col(stat="identity", position="dodge")
                + scale_fill_manual(values=("#79706e", "#86bcb6", "#F28e2b"))
                + theme(legend_title=element_blank())
            )
        else:
            p = (
                ggplot(res, aes(x="c_call", y="counts", fill="group"))
                + coord_flip()
                + theme_classic()
                + xlab("c_call")
                + ylab("% c calls")
                + geom_col(stat="identity", position="dodge")
                + scale_fill_manual(values=("#79706e", "#86bcb6"))
                + theme(legend_title=element_blank())
            )
        if save_plot:
            _file3 = "{}/assign_isotype.pdf".format(os.path.dirname(filePath))
            save_as_pdf_pages([p], filename=_file3)
            if show_plot:
                print(p)
        else:
            if show_plot:
                print(p)
    # move and rename
    move_to_tmp(fasta, filename_prefix)
    make_all(fasta, filename_prefix)
    rename_dandelion(fasta, filename_prefix)


def assign_isotypes(
    fastas: Sequence,
    fileformat: Literal["blast", "changeo", "airr"] = "blast",
    org: Literal["human", "mouse"] = "human",
    correct_c_call: bool = True,
    correction_dict: Optional[Dict] = None,
    plot: bool = True,
    save_plot: bool = False,
    show_plot: bool = True,
    figsize: Tuple[Union[int, float], Union[int, float]] = (4, 4),
    blastdb: Optional[str] = None,
    allele: bool = False,
    filename_prefix: Optional[Union[Sequence, str]] = None,
    verbose: bool = False,
):
    """
    Annotate contigs with constant region call using blastn.

    Parameters
    ----------
    fastas : Sequence
        list or sequence of paths to fasta files.
    fileformat : str
        format of V(D)J file/objects. Default is 'blast'. Also accepts 'changeo' (same behaviour as 'blast') and 'airr'.
    org : str
        organism of reference folder. Default is 'human'.
    correct_c_call : bool
        whether or not to adjust the c_calls after blast based on provided primers specified in `primer_dict` option.
        Default is True.
    correction_dict : Dict, Optional
        a nested dictionary contain isotype/c_genes as keys and primer sequences as records to use for correcting
        annotated c_calls. Defaults to a curated dictionary for human sequences if left as none.
    plot : bool
        whether or not to plot reassignment summary metrics. Default is True.
    save_plot : bool
        whether or not to save plots.
    show_plot : bool
        whether or not to show plots.
    figsize : Tuple[Union[int,float], Union[int,float]]
        size of figure. Default is (4, 4).
    blastdb : str, Optional
        path to blast database. Defaults to `$BLASTDB` environmental variable.
    allele : bool
        whether or not to return allele calls. Default is False.
    filename_prefix : str, Optional
        list of prefixes of file names preceding '_contig'. None defaults to 'filtered'.
    verbose : bool
        whether or not to print the blast command in terminal. Default is False.

    Returns
    -------
    V(D)J tsv files with constant genes annotated.
    """
    if type(fastas) is not list:
        fastas = [fastas]
    if type(filename_prefix) is not list:
        filename_prefix = [filename_prefix]
    if all(t is None for t in filename_prefix):
        filename_prefix = [None for f in fastas]

    if verbose:
        print("Assign isotypes \n")

    for i in range(0, len(fastas)):
        assign_isotype(
            fastas[i],
            fileformat=fileformat,
            org=org,
            correct_c_call=correct_c_call,
            correction_dict=correction_dict,
            plot=plot,
            save_plot=save_plot,
            show_plot=show_plot,
            figsize=figsize,
            blastdb=blastdb,
            allele=allele,
            filename_prefix=filename_prefix[i],
            verbose=verbose,
        )


def reannotate_genes(
    data: Sequence,
    igblast_db: Optional[str] = None,
    germline: Optional[Union[str, PathLike]] = None,
    org: Literal["human", "mouse"] = "human",
    loci: Literal["ig", "tr"] = "ig",
    extended: bool = True,
    filename_prefix: Optional[Union[Sequence, str]] = None,
    flavour: Literal["strict", "original"] = "strict",
    min_j_match: int = 7,
    min_d_match: int = 9,
    v_evalue: float = 1e-4,
    d_evalue: float = 1e-3,
    j_evalue: float = 1e-4,
    reassign_dj: bool = True,
    overwrite: bool = True,
    dust: Optional[Union[Literal["yes", "no"], str]] = "no",
    verbose: bool = False,
):
    """
    Reannotate cellranger fasta files with igblastn and parses to airr format.

    Parameters
    ----------
    data : Sequence
        list of fasta file locations, or folder name containing fasta files.
        if provided as a single string, it will first be converted to a list;
        this allows for the function to be run on single/multiple samples.
    igblast_db : str, PathLike, Optional
        path to igblast database folder. Defaults to `$IGDATA` environmental
        variable.
    germline : str, PathLike, Optional
        path to germline database folder. Defaults to `$GERMLINE` environmental
        variable.
    org : str
        organism of germline database. Default is 'human'.
    loci : str
        mode for igblastn. Default is 'ig' for BCRs. Also accepts 'tr' for
        TCRs.
    extended : bool
        whether or not to transfer additional 10X annotions to output file.
        Default is True.
    filename_prefix : str, Optional
        list of prefixes of file names preceding '_contig'. None defaults
        to 'filtered'.
    flavour : str
        Either 'dandelion' or 'immcantation'. Determines how igblastnshould
        be run. Running in 'dandelion' flavour will add the additional the
        evalue and min_d_match options to the run.
    v_evalue : float
        This is the statistical significance threshold for reporting matches
        against database sequences. Lower EXPECT thresholds are more stringent
        and report only high similarity matches. Choose higher EXPECT value
        (for example 1 or more) if you expect a low identity between your query
        sequence and the targets. for v gene.
    d_evalue : float
        This is the statistical significance threshold for reporting matches
        against database sequences. Lower EXPECT thresholds are more stringent
        and report only high similarity matches. Choose higher EXPECT value
        (for example 1 or more) if you expect a low identity between your query
        sequence and the targets. for d gene.
    j_evalue : float
        This is the statistical significance threshold for reporting matches
        against database sequences. Lower EXPECT thresholds are more stringent
        and report only high similarity matches. Choose higher EXPECT value
        (for example 1 or more) if you expect a low identity between your query
        sequence and the targets. for j gene.
    min_d_match : int
        Minimum D gene nucleotide matches. This controls the threshold for
        D gene detection. You can set the minimal number of required
        consecutive nucleotide matches between the query sequence and the D
        genes based on your own criteria. Note that the matches do not include
        overlapping matches at V-D or D-J junctions.
    min_j_match : int
        Minimum D gene nucleotide matches. This controls the threshold for
        D gene detection. You can set the minimal number of required
        consecutive nucleotide matches between the query sequence and the D
        genes based on your own criteria. Note that the matches do not include
        overlapping matches at V-D or D-J junctions.
    reassign_dj : bool
        whether or not to perform a targetted blastn reassignment for D and J genes.
        Default is False.
    dust: str
        dustmasker options. Filter query sequence with DUST
        Format: 'yes', or 'no' to disable. Accepts str.
        If None, defaults to `20 64 1`.
    overwrite: bool
        whether or not to overwrite the assignment if flavour = 'strict'.
    verbose :
        whether or not to print the igblast command used in the terminal.
        Default is False.

    Returns
    -------
    V(D)J data file in airr/changeo data format.
    """
    if type(data) is not list:
        data = [data]
    if type(filename_prefix) is not list:
        filename_prefix = [filename_prefix]
    if all(t is None for t in filename_prefix):
        filename_prefix = [None for d in data]

    filePath = None
    for i in tqdm(range(0, len(data)), desc="Assigning genes "):
        filePath = check_filepath(
            data[i], filename_prefix=filename_prefix[i], endswith=".fasta"
        )
        if filePath is None:
            if filename_prefix[i] is not None:
                raise FileNotFoundError(
                    "Path to fasta file with filename prefix `{}_contig` is unknown. ".format(
                        filename_prefix[i]
                    )
                    + "Please specify path to fasta file or folder containing fasta file."
                )
            else:
                raise FileNotFoundError(
                    "Path to fasta file is unknown. "
                    + "Please specify path to fasta file or folder containing fasta file."
                )

        if verbose:
            print("Processing {} \n".format(filePath))

        if flavour == "original":
            assigngenes_igblast(
                filePath,
                igblast_db=igblast_db,
                org=org,
                loci=loci,
                verbose=verbose,
            )
        elif flavour == "strict":
            run_igblastn(
                filePath,
                igblast_db=igblast_db,
                org=org,
                loci=loci,
                evalue=v_evalue,
                min_d_match=min_d_match,
                verbose=verbose,
            )
        makedb_igblast(
            filePath,
            org=org,
            germline=germline,
            extended=extended,
            verbose=verbose,
        )
        # block this for now, until I figure out if it's
        # worth it
        if flavour == "strict":
            if reassign_dj:
                assign_DJ(
                    fasta=filePath,
                    org=org,
                    loci=loci,
                    call="j",
                    database=igblast_db,
                    evalue=j_evalue,
                    filename_prefix=filename_prefix,
                    dust=dust,
                    word_size=min_j_match,
                    overwrite=overwrite,
                    verbose=verbose,
                )
                assign_DJ(
                    filePath,
                    org=org,
                    loci=loci,
                    call="d",
                    database=igblast_db,
                    evalue=d_evalue,
                    filename_prefix=filename_prefix,
                    dust=dust,
                    word_size=min_d_match,
                    overwrite=overwrite,
                    verbose=verbose,
                )

    if loci == "tr":
        change_file_location(data, filename_prefix)
        if flavour == "strict":
            mask_dj(data, filename_prefix, d_evalue, j_evalue)
        move_to_tmp(data, filename_prefix)
        make_all(data, filename_prefix)
        rename_dandelion(data, filename_prefix, endswith="_igblast_db-pass.tsv")


def reassign_alleles(
    data: Sequence,
    combined_folder: Union[str, PathLike],
    v_germline: Optional[str] = None,
    germline: Optional[Union[str, PathLike]] = None,
    org: Literal["human", "mouse"] = "human",
    v_field: Literal["v_call", "v_call_genotyped"] = "v_call_genotyped",
    germ_types: Literal["full", "dmask", "vonly", "regions"] = "dmask",
    novel: bool = True,
    cloned: bool = False,
    plot: bool = True,
    save_plot: bool = False,
    show_plot: bool = True,
    figsize: Tuple[Union[int, float], Union[int, float]] = (4, 3),
    sample_id_dictionary: Optional[Dict] = None,
    filename_prefix: Optional[Union[Sequence, str]] = None,
    verbose: bool = False,
):
    """
    Correct allele calls based on a personalized genotype using tigger.

    It uses a subject-specific genotype to correct correct preliminary allele
    assignments of a set ofsequences derived from a single subject.

    Parameters
    ----------
    data : Sequence
        list of data folders containing the .tsv files. if provided as a single
        string, it will first be converted to a list; this allows for the
        function to be run on single/multiple samples.
    combined_folder : str, PathLike
        name of folder for concatenated data file and genotyped files.
    v_germline : str, Optional
        path to heavy chain v germline fasta. Defaults to IGHV fasta in
        `$GERMLINE` environmental variable.
    germline : str, Optional
        path to germline database folder. Defaults to `$GERMLINE` environmental
        variable.
    org : str
        organism of germline database. Default is 'human'.
    v_field : str
        name of column containing the germline V segment call. Default is
        v_call_genotyped' (airr) for after tigger.
    germ_types : str
        Specify type of germline for reconstruction. Accepts one of : 'full',
        'dmask', 'vonly', 'region'. Default is 'dmask'.
    novel : bool
        whether or not to run novel allele discovery during tigger-genotyping.
        Default is True (yes).
    cloned : bool
        whether or not to run CreateGermlines.py with `--cloned`.
    plot : bool
        whether or not to plot reassignment summary metrics. Default is True.
    save_plot : bool
        whether or not to save plot.
    show_plot : bool
        whether or not to show plot.
    figsize : Tuple[Union[int,float], Union[int,float]]
        size of figure. Default is (4, 3).
    sample_id_dictionary : dict, Optional
        dictionary for creating a sample_id column in the concatenated file.
    filename_prefix : str, Optional
        list of prefixes of file names preceding '_contig'. None defaults to
        'filtered'.
    verbose : bool
        Whether or not to print the command used in the terminal. Default is
        False.

    Returns
    -------
    Individual V(D)J data files with v_call_genotyped column containing
    reassigned heavy chain v calls
    """
    fileformat = "blast"
    if type(data) is not list:
        data = [data]
    if type(filename_prefix) is not list:
        filename_prefix = [filename_prefix]
    if all(t is None for t in filename_prefix):
        filename_prefix = [None for d in data]

    informat_dict = {
        "changeo": "_igblast_db-pass.tsv",
        "blast": "_igblast_db-pass.tsv",
        "airr": "_igblast_gap.tsv",
    }
    germpass_dict = {
        "changeo": "_igblast_db-pass_germ-pass.tsv",
        "blast": "_igblast_db-pass_germ-pass.tsv",
        "airr": "_igblast_gap_germ-pass.tsv",
    }
    heavy_dict = {
        "changeo": "_igblast_db-pass_heavy_parse-select.tsv",
        "blast": "_igblast_db-pass_heavy_parse-select.tsv",
        "airr": "_igblast_gap_heavy_parse-select.tsv",
    }
    light_dict = {
        "changeo": "_igblast_db-pass_light_parse-select.tsv",
        "blast": "_igblast_db-pass_light_parse-select.tsv",
        "airr": "_igblast_gap_light_parse-select.tsv",
    }
    fileformat_dict = {
        "changeo": "_igblast_db-pass_genotyped.tsv",
        "blast": "_igblast_db-pass_genotyped.tsv",
        "airr": "_igblast_gap_genotyped.tsv",
    }
    fileformat_passed_dict = {
        "changeo": "_igblast_db-pass_genotyped_germ-pass.tsv",
        "blast": "_igblast_db-pass_genotyped_germ-pass.tsv",
        "airr": "_igblast_gap_genotyped_germ-pass.tsv",
    }
    inferred_fileformat_dict = {
        "changeo": "_igblast_db-pass_inferredGenotype.txt",
        "blast": "_igblast_db-pass_inferredGenotype.txt",
        "airr": "_igblast_gap_inferredGenotype.txt",
    }
    germline_dict = {
        "changeo": "_igblast_db-pass_genotype.fasta",
        "blast": "_igblast_db-pass_genotype.fasta",
        "airr": "_igblast_gap_genotype.fasta",
    }
    fform_dict = {"blast": "airr", "airr": "airr", "changeo": "changeo"}

    filepathlist_heavy = []
    filepathlist_light = []
    filePath = None
    sampleNames_dict = {}
    filePath_dict = {}
    for i in tqdm(range(0, len(data)), desc="Processing data file(s) "):
        filePath = check_filepath(
            data[i],
            filename_prefix=filename_prefix[i],
            endswith=informat_dict[fileformat],
            subdir="tmp",
        )
        if filePath is None:
            raise FileNotFoundError(
                "Path to .tsv file for {} is unknown. ".format(data[i])
                + "Please specify path to reannotated .tsv file or folder "
                + "containing reannotated .tsv file."
            )

        filePath_heavy = filePath.replace(
            informat_dict[fileformat], heavy_dict[fileformat]
        )
        filePath_light = filePath.replace(
            informat_dict[fileformat], light_dict[fileformat]
        )

        if sample_id_dictionary is not None:
            sampleNames_dict[filePath] = sample_id_dictionary[data[i]]
        else:
            sampleNames_dict[filePath] = str(data[i])

        filePath_dict[str(data[i])] = filePath

        # splitting up to heavy chain and light chain files
        parsedb_heavy(filePath)
        parsedb_light(filePath)

        # add to counter
        filepathlist_heavy.append(filePath_heavy)
        filepathlist_light.append(filePath_light)

    # make output directory
    outDir = combined_folder.rstrip("/")
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    # concatenate
    if len(filepathlist_heavy) > 1:
        print("Concatenating objects")
        cmd1 = " ".join(
            [
                'awk "FNR==1 && NR!=1 { while (/^sequence_id/) getline; } 1 {print}"'
            ]
            + [f for f in filepathlist_heavy]
            + [">"]
            + [outDir + "/" + outDir + "_heavy" + informat_dict[fileformat]]
        )
        cmd2 = " ".join(
            [
                'awk "FNR==1 && NR!=1 { while (/^sequence_id/) getline; } 1 {print}"'
            ]
            + [f for f in filepathlist_light]
            + [">"]
            + [outDir + "/" + outDir + "_light" + informat_dict[fileformat]]
        )
    else:
        cmd1 = " ".join(
            [
                'awk "FNR==1 && NR!=1 { while (/^sequence_id/) getline; } 1 {print}"'
            ]
            + [filepathlist_heavy[0]]
            + [">"]
            + [outDir + "/" + outDir + "_heavy" + informat_dict[fileformat]]
        )
        cmd2 = " ".join(
            [
                'awk "FNR==1 && NR!=1 { while (/^sequence_id/) getline; } 1 {print}"'
            ]
            + [filepathlist_light[0]]
            + [">"]
            + [outDir + "/" + outDir + "_light" + informat_dict[fileformat]]
        )

    if verbose:
        print("Running command: %s\n" % (cmd1))
        print("Running command: %s\n" % (cmd2))
    os.system(cmd1)
    os.system(cmd2)

    novel_dict = {True: "YES", False: "NO"}
    if novel:
        try:
            print("      Running tigger-genotype with novel allele discovery.")
            tigger_genotype(
                outDir + "/" + outDir + "_heavy" + informat_dict[fileformat],
                v_germline=v_germline,
                fileformat=fform_dict[fileformat],
                novel_=novel_dict[novel],
                verbose=verbose,
            )
            creategermlines(
                outDir + "/" + outDir + "_heavy" + fileformat_dict[fileformat],
                germtypes=germ_types,
                mode="heavy",
                genotype_fasta=outDir
                + "/"
                + outDir
                + "_heavy"
                + germline_dict[fileformat],
                germline=germline,
                v_field=v_field,
                verbose=verbose,
                cloned=cloned,
            )
            _ = load_data(
                outDir
                + "/"
                + outDir
                + "_heavy"
                + fileformat_passed_dict[fileformat]
            )
        except:
            try:
                print("      Novel allele discovery execution halted.")
                print(
                    "      Attempting to run tigger-genotype without novel allele discovery."
                )
                tigger_genotype(
                    outDir
                    + "/"
                    + outDir
                    + "_heavy"
                    + informat_dict[fileformat],
                    v_germline=v_germline,
                    fileformat=fform_dict[fileformat],
                    novel_=novel_dict[False],
                    verbose=verbose,
                )
                creategermlines(
                    outDir
                    + "/"
                    + outDir
                    + "_heavy"
                    + fileformat_dict[fileformat],
                    germtypes=germ_types,
                    mode="heavy",
                    genotype_fasta=outDir
                    + "/"
                    + outDir
                    + "_heavy"
                    + germline_dict[fileformat],
                    germline=germline,
                    v_field=v_field,
                    verbose=verbose,
                    cloned=cloned,
                )
                _ = load_data(
                    outDir
                    + "/"
                    + outDir
                    + "_heavy"
                    + fileformat_passed_dict[fileformat]
                )
            except:
                print(
                    "     Insufficient contigs for running tigger-genotype. Defaulting to original heavy chain v_calls."
                )
                tigger_failed = ""
    else:
        try:
            print(
                "      Running tigger-genotype without novel allele discovery."
            )
            tigger_genotype(
                outDir + "/" + outDir + "_heavy" + informat_dict[fileformat],
                v_germline=v_germline,
                fileformat=fform_dict[fileformat],
                novel_=novel_dict[False],
                verbose=verbose,
            )
            creategermlines(
                outDir + "/" + outDir + "_heavy" + fileformat_dict[fileformat],
                germtypes=germ_types,
                mode="heavy",
                genotype_fasta=outDir
                + "/"
                + outDir
                + "_heavy"
                + germline_dict[fileformat],
                germline=germline,
                v_field=v_field,
                verbose=verbose,
                cloned=cloned,
            )
            _ = load_data(
                outDir
                + "/"
                + outDir
                + "_heavy"
                + fileformat_passed_dict[fileformat]
            )
        except:
            print(
                "      Insufficient contigs for running tigger-genotype. Defaulting to original heavy chain v_calls."
            )
            tigger_failed = ""

    if "tigger_failed" in locals():
        creategermlines(
            outDir + "/" + outDir + "_heavy" + informat_dict[fileformat],
            germtypes=germ_types,
            mode="heavy",
            genotype_fasta=None,
            germline=germline,
            v_field="v_call",
            verbose=verbose,
            cloned=cloned,
        )
        creategermlines(
            outDir + "/" + outDir + "_light" + informat_dict[fileformat],
            germtypes=germ_types,
            mode="light",
            genotype_fasta=None,
            germline=germline,
            v_field="v_call",
            verbose=verbose,
            cloned=cloned,
        )
        print(
            "      For convenience, entries for heavy chain in `v_call` are copied to `v_call_genotyped`."
        )
        heavy = load_data(
            outDir + "/" + outDir + "_heavy" + germpass_dict[fileformat]
        )
        heavy["v_call_genotyped"] = heavy["v_call"]
        print(
            "      For convenience, entries for light chain `v_call` are copied to `v_call_genotyped`."
        )
        light = load_data(
            outDir + "/" + outDir + "_light" + germpass_dict[fileformat]
        )
        light["v_call_genotyped"] = light["v_call"]
    else:
        creategermlines(
            outDir + "/" + outDir + "_light" + informat_dict[fileformat],
            germtypes=germ_types,
            mode="light",
            genotype_fasta=None,
            germline=germline,
            v_field="v_call",
            verbose=verbose,
            cloned=cloned,
        )
        heavy = load_data(
            outDir
            + "/"
            + outDir
            + "_heavy"
            + fileformat_passed_dict[fileformat]
        )
        print(
            "      For convenience, entries for light chain `v_call` are copied to `v_call_genotyped`."
        )
        light = load_data(
            outDir + "/" + outDir + "_light" + germpass_dict[fileformat]
        )
        light["v_call_genotyped"] = light["v_call"]

    sampledict = {}
    heavy["sample_id"], light["sample_id"] = None, None
    for file in sampleNames_dict.keys():
        dat_f = load_data(file)
        dat_f["sample_id"] = sampleNames_dict[file]
        heavy["sample_id"].update(dat_f["sample_id"])
        light["sample_id"].update(dat_f["sample_id"])

    dat_ = heavy.append(light)
    if "cell_id" in dat_.columns:
        dat_.sort_values(by="cell_id", inplace=True)
    else:
        dat_.sort_values(by="sequence_id", inplace=True)

    if plot:
        if "tigger_failed" not in locals():
            print("Returning summary plot")
            inferred_genotype = (
                outDir
                + "/"
                + outDir
                + "_heavy"
                + inferred_fileformat_dict[fileformat]
            )
            inf_geno = pd.read_csv(inferred_genotype, sep="\t", dtype="object")

            s2 = set(inf_geno["gene"])
            results = []
            try:
                for samp in list(set(heavy["sample_id"])):
                    res_x = heavy[(heavy["sample_id"] == samp)]
                    V_ = [
                        re.sub("[*][0-9][0-9]", "", v) for v in res_x["v_call"]
                    ]
                    V_g = [
                        re.sub("[*][0-9][0-9]", "", v)
                        for v in res_x["v_call_genotyped"]
                    ]
                    s1 = set(
                        list(
                            ",".join(
                                [",".join(list(set(v.split(",")))) for v in V_]
                            ).split(",")
                        )
                    )
                    setdiff = s1 - s2
                    ambiguous = (
                        ["," in i for i in V_].count(True) / len(V_) * 100,
                        ["," in i for i in V_g].count(True) / len(V_g) * 100,
                    )
                    not_in_genotype = (
                        [i in setdiff for i in V_].count(True) / len(V_) * 100,
                        [i in setdiff for i in V_g].count(True)
                        / len(V_g)
                        * 100,
                    )
                    stats = pd.DataFrame(
                        [ambiguous, not_in_genotype],
                        columns=["ambiguous", "not_in_genotype"],
                        index=["before", "after"],
                    ).T
                    stats.index.set_names(["vgroup"], inplace=True)
                    stats.reset_index(drop=False, inplace=True)
                    stats["sample_id"] = samp
                    # stats['donor'] = str(combined_folder)
                    results.append(stats)
                results = pd.concat(results)
                ambiguous_table = results[results["vgroup"] == "ambiguous"]
                not_in_genotype_table = results[
                    results["vgroup"] == "not_in_genotype"
                ]
                ambiguous_table.reset_index(inplace=True, drop=True)
                not_in_genotype_table.reset_index(inplace=True, drop=True)
                # melting the dataframe
                ambiguous_table_before = ambiguous_table.drop("after", axis=1)
                ambiguous_table_before.rename(
                    columns={"before": "var"}, inplace=True
                )
                ambiguous_table_before["var_group"] = "before"
                ambiguous_table_after = ambiguous_table.drop("before", axis=1)
                ambiguous_table_after.rename(
                    columns={"after": "var"}, inplace=True
                )
                ambiguous_table_after["var_group"] = "after"
                ambiguous_table = pd.concat(
                    [ambiguous_table_before, ambiguous_table_after]
                )
                not_in_genotype_table_before = not_in_genotype_table.drop(
                    "after", axis=1
                )
                not_in_genotype_table_before.rename(
                    columns={"before": "var"}, inplace=True
                )
                not_in_genotype_table_before["var_group"] = "before"
                not_in_genotype_table_after = not_in_genotype_table.drop(
                    "before", axis=1
                )
                not_in_genotype_table_after.rename(
                    columns={"after": "var"}, inplace=True
                )
                not_in_genotype_table_after["var_group"] = "after"
                not_in_genotype_table = pd.concat(
                    [not_in_genotype_table_before, not_in_genotype_table_after]
                )
                ambiguous_table["var_group"] = ambiguous_table[
                    "var_group"
                ].astype("category")
                not_in_genotype_table["var_group"] = not_in_genotype_table[
                    "var_group"
                ].astype("category")
                ambiguous_table["var_group"].cat.reorder_categories(
                    ["before", "after"], inplace=True
                )
                not_in_genotype_table["var_group"].cat.reorder_categories(
                    ["before", "after"], inplace=True
                )

                options.figure_size = figsize
                final_table = pd.concat(
                    [ambiguous_table, not_in_genotype_table]
                )
                p = (
                    ggplot(
                        final_table,
                        aes(x="sample_id", y="var", fill="var_group"),
                    )
                    + coord_flip()
                    + theme_classic()
                    + xlab("sample_id")
                    + ylab("% allele calls")
                    + ggtitle("Genotype reassignment with TIgGER")
                    + geom_bar(stat="identity")
                    + facet_grid("~" + str("vgroup"), scales="free_y")
                    + scale_fill_manual(values=("#86bcb6", "#F28e2b"))
                    + theme(legend_title=element_blank())
                )
                if save_plot:
                    savefile = outDir + "/" + outDir + "_reassign_alleles.pdf"
                    save_as_pdf_pages([p], filename=savefile)
                    if show_plot:
                        print(p)
                else:
                    if show_plot:
                        print(p)
            except:
                print("Error in plotting encountered. Skipping.")
                pass
        else:
            pass
    sleep(0.5)
    # if split_write_out:
    if "tigger_failed" in locals():
        print(
            "Although tigger-genotype was not run successfully, file will still be saved with `_genotyped.tsv`"
            "extension for convenience."
        )
    for s in tqdm(data, desc="Writing out to individual folders "):
        if sample_id_dictionary is not None:
            out_file = dat_[dat_["sample_id"] == sample_id_dictionary[s]]
        else:
            out_file = dat_[dat_["sample_id"] == s]
        outfilepath = filePath_dict[s]
        write_airr(out_file, outfilepath.replace(".tsv", "_genotyped.tsv"))


def create_germlines(
    self: Union[Dandelion, pd.DataFrame, str, PathLike],
    germline: Optional[Union[str, PathLike]] = None,
    org: Literal["human", "mouse"] = "human",
    seq_field: Literal["sequence_alignment"] = "sequence_alignment",
    v_field: Literal["v_call", "v_call_genotyped"] = "v_call",
    d_field: Literal["d_call"] = "d_call",
    j_field: Literal["j_call"] = "j_call",
    germ_types: Literal["full", "dmask", "vonly", "regions"] = "dmask",
    fileformat: Literal["changeo", "airr"] = "airr",
    initialize_metadata: bool = False,
) -> Dandelion:
    """
    Run CreateGermlines.py to reconstruct the germline V(D)J sequence.

    Run CreateGermlines.py to reconstruct the germline V(D)J sequence.

    Parameters
    ----------
    self : Dandelion, pd.DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr
        file after clones have been determined.
    germline : str, Optional
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
        Specify type(s) of germlines to include full germline, germline with D segment masked,
        or germline for V segment only. Default is 'dmask'.
    fileformat : str
        format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.

    Returns
    -------
    V(D)J data file with reconstructed germline sequences.
    """
    start = logg.info("Reconstructing germline sequences")
    env = os.environ.copy()

    if self.__class__ != Dandelion:
        if germline is None:
            try:
                gml = env["GERMLINE"]
            except:
                raise KeyError(
                    "Environmental variable GERMLINE must be set. "
                    + "Otherwise, please provide path to folder containing germline fasta files."
                )
            gml = gml + "imgt/" + org + "/vdj/"
        else:
            if os.path.isdir(germline):
                env["GERMLINE"] = germline
                gml = germline
            elif type(germline) is list:
                gml = germline
            elif type(germline) is dict:
                gml = germline
    else:
        if len(self.germline) == 0:
            if germline is None:
                try:
                    gml = env["GERMLINE"]
                except:
                    raise KeyError(
                        "Environmental variable GERMLINE must be set. "
                        + "Otherwise, please provide path to folder containing germline fasta files."
                    )
                gml = gml + "imgt/" + org + "/vdj/"
            else:
                if os.path.isdir(germline):
                    env["GERMLINE"] = germline
                    gml = germline
                elif type(germline) is list:
                    gml = germline
                elif type(germline) is dict:
                    gml = germline

    def _parseChangeO(record):
        """
        Parse a dictionary to a Receptor object.

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
        Parse a dictionary of AIRR records to a Receptor object.

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

    def _create_germlines_object(
        self,
        references,
        seq_field,
        v_field,
        d_field,
        j_field,
        germ_types,
        fileformat,
    ):
        """
        Write germline sequences to tab-delimited database file.

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
            raise ValueError("Invalid format %s" % fileformat)

        # Define output germline fields
        germline_fields = OrderedDict()
        seq_type = seq_field.split("_")[-1]
        if "full" in germ_types:
            germline_fields["full"] = "germline_" + seq_type
        if "dmask" in germ_types:
            germline_fields["dmask"] = "germline_" + seq_type + "_d_mask"
        if "vonly" in germ_types:
            germline_fields["vonly"] = "germline_" + seq_type + "_v_region"
        if "regions" in germ_types:
            germline_fields["regions"] = "germline_regions"

        if type(references) is dict:
            reference_dict = references
        else:
            if type(references) is not list:
                ref = [references]
            else:
                ref = references
            reference_dict = readGermlines(ref)
        # Check for IMGT-gaps in germlines
        if all("..." not in x for x in reference_dict.values()):
            warnings.warn(
                UserWarning(
                    "Germline reference sequences do not appear to contain IMGT-numbering spacers."
                    + " Results may be incorrect."
                )
            )

        required = [
            "v_germ_start_imgt",
            "d_germ_start",
            "j_germ_start",
            "np1_length",
            "np2_length",
        ]

        if self.__class__ == Dandelion:
            if isinstance(self.data, pd.DataFrame):
                # Check for required columns
                try:
                    checkFields(required, self.data.columns, schema=schema)
                except LookupError as e:
                    print(e)

                # Count input
                # total_count = len(self.data)

                # Check for existence of fields
                for f in [v_field, d_field, j_field, seq_field]:
                    if f not in self.data.columns:
                        raise NameError(
                            "%s field does not exist in input database file."
                            % f
                        )
                # Translate to Receptor attribute names
                v_field_ = schema.toReceptor(v_field)
                d_field_ = schema.toReceptor(d_field)
                j_field_ = schema.toReceptor(j_field)
                seq_field_ = schema.toReceptor(seq_field)
                # clone_field = schema.toReceptor(clone_field)

                # Define Receptor iterator
                receptor_iter = (
                    (
                        self.data.loc[
                            x,
                        ].sequence_id,
                        self.data.loc[
                            x,
                        ],
                    )
                    for x in self.data.index
                )

            else:
                raise LookupError(
                    "Please initialise the Dandelion object with a dataframe in data slot."
                )
        elif self.__class__ == pd.DataFrame:
            try:
                checkFields(required, self.columns, schema=schema)
            except LookupError as e:
                print(e)

            # Count input
            # total_count = len(self)
            # Check for existence of fields
            for f in [v_field, d_field, j_field, seq_field]:
                if f not in self.columns:
                    raise NameError(
                        "%s field does not exist in input database file." % f
                    )
            # Translate to Receptor attribute names
            v_field_ = schema.toReceptor(v_field)
            d_field_ = schema.toReceptor(d_field)
            j_field_ = schema.toReceptor(j_field)
            seq_field_ = schema.toReceptor(seq_field)
            # clone_field = schema.toReceptor(clone_field)
            # Define Receptor iterator
            receptor_iter = (
                (
                    self.loc[
                        x,
                    ].sequence_id,
                    self.loc[
                        x,
                    ],
                )
                for x in self.index
            )

        out = {}
        # Iterate over rows
        for key, records in tqdm(
            receptor_iter,
            desc="   Building {} germline sequences".format(germ_types),
        ):
            # Define iteration variables
            # Build germline for records
            if fileformat == "airr":
                germ_log, glines, genes = buildGermline(
                    _parseAIRR(dict(records)),
                    reference_dict,
                    seq_field=seq_field_,
                    v_field=v_field_,
                    d_field=d_field_,
                    j_field=j_field_,
                )
            elif fileformat == "changeo":
                germ_log, glines, genes = buildGermline(
                    _parseChangeO(dict(records)),
                    reference_dict,
                    seq_field=seq_field_,
                    v_field=v_field_,
                    d_field=d_field_,
                    j_field=j_field_,
                )
            else:
                raise AttributeError(
                    "%s is not acceptable file format." % fileformat
                )

            if glines is not None:
                # Add glines to Receptor record
                annotations = {}
                if "full" in germ_types:
                    annotations[germline_fields["full"]] = glines["full"]
                if "dmask" in germ_types:
                    annotations[germline_fields["dmask"]] = glines["dmask"]
                if "vonly" in germ_types:
                    annotations[germline_fields["vonly"]] = glines["vonly"]
                if "regions" in germ_types:
                    annotations[germline_fields["regions"]] = glines["regions"]
                out.update({key: annotations})
        germline_df = pd.DataFrame.from_dict(out, orient="index")

        if self.__class__ == Dandelion:
            # datx = load_data(self.data)
            for x in germline_df.columns:
                self.data[x] = pd.Series(germline_df[x])

        elif self.__class__ == pd.DataFrame:
            datx = load_data(self)
            for x in germline_df.columns:
                datx[x] = pd.Series(germline_df[x])
            try:
                output = Dandelion(
                    data=datx, germline=reference_dict, initialize=True
                )
            except:
                output = Dandelion(
                    data=datx, germline=reference_dict, initialize=False
                )
            return output
        sleep(0.5)
        logg.info(
            " finished",
            time=start,
            deep=(
                "Updated Dandelion object: \n"
                "   'data', updated germline alignment in contig-indexed clone table\n"
                "   'germline', updated germline reference\n"
            ),
        )

    def _create_germlines_file(
        file,
        references,
        seq_field,
        v_field,
        d_field,
        j_field,
        germ_types,
        fileformat,
    ):
        """
        Write germline sequences to tab-delimited database file.

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
            raise ValueError("Invalid format %s" % fileformat)

        # Define output germline fields
        germline_fields = OrderedDict()
        seq_type = seq_field.split("_")[-1]
        if "full" in germ_types:
            germline_fields["full"] = "germline_" + seq_type
        if "dmask" in germ_types:
            germline_fields["dmask"] = "germline_" + seq_type + "_d_mask"
        if "vonly" in germ_types:
            germline_fields["vonly"] = "germline_" + seq_type + "_v_region"
        if "regions" in germ_types:
            germline_fields["regions"] = "germline_regions"

        if type(references) is dict:
            reference_dict = references
        else:
            if type(references) is not list:
                ref = [references]
            else:
                ref = references
            reference_dict = readGermlines(ref)
        # Check for IMGT-gaps in germlines
        if all("..." not in x for x in reference_dict.values()):
            warnings.warn(
                UserWarning(
                    "Germline reference sequences do not appear to contain IMGT-numbering spacers. "
                    + "Results may be incorrect."
                )
            )

        required = [
            "v_germ_start_imgt",
            "d_germ_start",
            "j_germ_start",
            "np1_length",
            "np2_length",
        ]

        # Get repertoire and open Db reader
        db_handle = open(file, "rt")
        db_iter = reader(db_handle)
        # Check for required columns
        try:
            checkFields(required, db_iter.fields, schema=schema)
        except LookupError as e:
            print(e)
        # Count input
        # total_count = countDbFile(file)
        # Check for existence of fields
        for f in [v_field, d_field, j_field, seq_field]:
            if f not in db_iter.fields:
                raise NameError(
                    "%s field does not exist in input database file." % f
                )
        # Translate to Receptor attribute names
        v_field_ = schema.toReceptor(v_field)
        d_field_ = schema.toReceptor(d_field)
        j_field_ = schema.toReceptor(j_field)
        seq_field_ = schema.toReceptor(seq_field)
        # clone_field = schema.toReceptor(clone_field)
        # Define Receptor iterator
        receptor_iter = ((x.sequence_id, [x]) for x in db_iter)

        out = {}
        # Iterate over rows
        for key, records in tqdm(
            receptor_iter,
            desc="   Building {} germline sequences".format(germ_types),
        ):
            # Define iteration variables
            # Build germline for records
            # if not isinstance(self.data, pd.DataFrame):
            records = list(records)
            germ_log, glines, genes = buildGermline(
                records[0],
                reference_dict,
                seq_field=seq_field_,
                v_field=v_field_,
                d_field=d_field_,
                j_field=j_field_,
            )
            if glines is not None:
                # Add glines to Receptor record
                annotations = {}
                if "full" in germ_types:
                    annotations[germline_fields["full"]] = glines["full"]
                if "dmask" in germ_types:
                    annotations[germline_fields["dmask"]] = glines["dmask"]
                if "vonly" in germ_types:
                    annotations[germline_fields["vonly"]] = glines["vonly"]
                if "regions" in germ_types:
                    annotations[germline_fields["regions"]] = glines["regions"]
                out.update({key: annotations})
        germline_df = pd.DataFrame.from_dict(out, orient="index")

        try:
            out = Dandelion(data=file, germline=reference_dict, initialize=True)
        except:
            out = Dandelion(
                data=file, germline=reference_dict, initialize=False
            )
        for x in germline_df.columns:
            out.data[x] = pd.Series(germline_df[x])

        if os.path.isfile(str(file)):
            out.write_airr(
                "{}/{}_germline_{}.tsv".format(
                    os.path.dirname(file),
                    os.path.basename(file).split(".tsv")[0],
                    germ_types,
                )
            )
        return out

    if (type(germline) is dict) or (type(germline) is list):
        if self.__class__ == Dandelion:
            _create_germlines_object(
                self,
                germline,
                seq_field,
                v_field,
                d_field,
                j_field,
                germ_types,
                fileformat,
            )
        elif self.__class__ == pd.DataFrame:
            return _create_germlines_object(
                self,
                germline,
                seq_field,
                v_field,
                d_field,
                j_field,
                germ_types,
                fileformat,
            )
        else:
            return _create_germlines_file(
                self,
                germline,
                seq_field,
                v_field,
                d_field,
                j_field,
                germ_types,
                fileformat,
            )
    else:
        if self.__class__ == Dandelion:
            if len(self.germline) != 0:
                _create_germlines_object(
                    self,
                    self.germline,
                    seq_field,
                    v_field,
                    d_field,
                    j_field,
                    germ_types,
                    fileformat,
                )
            else:
                _create_germlines_object(
                    self,
                    gml,
                    seq_field,
                    v_field,
                    d_field,
                    j_field,
                    germ_types,
                    fileformat,
                )
        elif self.__class__ == pd.DataFrame:
            return _create_germlines_object(
                self,
                gml,
                seq_field,
                v_field,
                d_field,
                j_field,
                germ_types,
                fileformat,
            )
        else:
            return _create_germlines_file(
                self,
                gml,
                seq_field,
                v_field,
                d_field,
                j_field,
                germ_types,
                fileformat,
            )


def filter_contigs(
    data: Union[Dandelion, pd.DataFrame, str],
    adata: Optional[AnnData] = None,
    filter_contig: bool = True,
    filter_rna: bool = False,
    filter_poorqualitycontig: bool = False,
    keep_highest_umi: bool = True,
    umi_foldchange_cutoff: int = 2,
    filter_vj_chains: bool = True,
    filter_missing: bool = True,
    productive_only: bool = True,
    simple: bool = False,
    save: Optional[str] = None,
    **kwargs,
) -> Tuple[Dandelion, AnnData]:
    """
    Filter doublets and poor quality cells and corresponding contigs based on provided V(D)J `DataFrame` and `AnnData`.

    Depends on a `AnnData`.obs slot populated with 'filter_rna' column. If the aligned sequence is an exact match
    between contigs, the contigs will be merged into the one with the highest umi count, adding the summing the
    umi count of the duplicated contigs to duplicate_count column. After this check, if there are still multiple
    contigs, cells with multiple contigs are filtered unless `keep_highest_umi` is False, where by the umi counts for
    each contig will then be compared and only the highest is retained. The contig with the highest umi that is
    > umi_foldchange_cutoff (default is empirically set at 2) will be retained. For productive heavy/long chains,
    if there are multiple contigs that survive the umi testing, then all contigs will be filtered. The default behaviour
    is to also filter cells with multiple light/short chains but this may sometimes be a true biological occurrence;
    toggling filter_vj_chains to False will rescue the mutltiplet light chains. Lastly, contigs with no corresponding
    cell barcode in the AnnData object is filtered if filter_missing is True. However, this may be useful to toggle to
    False if more contigs are preferred to be kept or for integrating with bulk reperotire seq data.

    Parameters
    ----------
    data : Dandeion, pd.DataDrame, str
        V(D)J airr/changeo data to filter. Can be pandas `DataFrame` object or file path as string.
    adata : AnnData, Optional
        AnnData object to filter. If not provided, will assume to keep all cells in the airr table.
    filter_contig : bool
        If True, V(D)J `DataFrame` object returned will be filtered. Default is True.
    filter_rna : bool
        If True, `AnnData` object returned will be filtered based on potential V(D)J doublets. Default is False.
    filter_poorqualitycontig : bool
        If True, barcodes marked with poor quality contigs will be filtered. Default is False; only relevant contigs are
        removed and RNA barcodes are kept.
    keep_highest_umi : bool
        If True, rescues IGH contigs with highest umi counts with a requirement that it passes the
        `umi_foldchange_cutoff` option. In addition, the sum of the all the heavy chain contigs must be greater than 3
        umi or all contigs will be filtered. Default is True.
    umi_foldchange_cutoff : int
        related to minimum fold change required to rescue heavy chain contigs/barcode otherwise they will be marked as
        doublets. Default is empirically set at 2-fold.
    filter_vj_chains : bool
        cells with multiple light chains will be marked to filter. Default is True.
    filter_missing : bool
        cells in V(D)J data not found in `AnnData` object will be marked to filter. Default is True. This may be useful
        for toggling to False if integrating with bulk data.
    productive_only : bool
        whether or not to retain only productive contigs.
    simple : bool
        simple filtering mode where only checks for potential gene assignment mismatches.
    save : str, Optional
        Only used if a pandas dataframe or dandelion object is provided. Specifying will save the formatted vdj table.
    **kwargs
        additional kwargs passed to `Dandelion.Dandelion`.

    Returns
    -------
    V(D)J `DataFrame` object in airr/changeo format and `AnnData` object.
    """
    start = logg.info("Filtering BCRs")
    if data.__class__ == Dandelion:
        dat_ = load_data(data.data)
    else:
        dat_ = load_data(data)

    if not simple:
        if productive_only:
            dat = dat_[dat_["productive"].isin(TRUES)].copy()
        else:
            dat = dat_.copy()
    else:
        dat = dat_.copy()

    if "cell_id" not in dat.columns:
        raise AttributeError(
            "VDJ data does not contain 'cell_id' column. Please make sure this is populated before filtering."
        )

    barcode = list(set(dat["cell_id"]))

    if adata is not None:
        adata_provided = True
        adata_ = adata.copy()
        if "filter_rna" not in adata_.obs:
            adata_.obs["filter_rna"] = "False"
        contig_check = pd.DataFrame(index=adata_.obs_names)
        bc_ = {}
        for b in barcode:
            bc_.update({b: "True"})
        contig_check["has_contig"] = pd.Series(bc_)
        contig_check.replace(np.nan, "No_contig", inplace=True)
        adata_.obs["has_contig"] = pd.Series(contig_check["has_contig"])
    else:
        adata_provided = False
        obs = pd.DataFrame(index=barcode)
        adata_ = ad.AnnData(obs=obs)
        adata_.obs["filter_rna"] = "False"
        adata_.obs["has_contig"] = "True"

    # rather than leaving a nan cell, i will create a 0 column for now
    if "duplicate_count" in dat and "umi_count" not in dat:
        dat["umi_count"] = dat["duplicate_count"]  # just do a simple swap?
    elif "duplicate_count" not in dat and "umi_count" in dat:
        dat["duplicate_count"] = dat["umi_count"]
    elif "duplicate_count" in dat and "umi_count" in dat:
        dat["umi_count"] = dat["duplicate_count"]

    if not simple:
        tofilter = FilterContigs(
            dat,
            keep_highest_umi,
            umi_foldchange_cutoff,
            filter_poorqualitycontig,
        )
    else:
        tofilter = FilterContigsLite(dat)

    poor_qual = tofilter.poor_qual.copy()
    h_doublet = tofilter.h_doublet.copy()
    l_doublet = tofilter.l_doublet.copy()
    drop_contig = tofilter.drop_contig.copy()
    umi_adjustment = tofilter.umi_adjustment.copy()

    if len(umi_adjustment) > 0:
        dat["duplicate_count"].update(umi_adjustment)

    poorqual = {c: "False" for c in adata_.obs_names}
    hdoublet = {c: "False" for c in adata_.obs_names}
    ldoublet = {c: "False" for c in adata_.obs_names}

    poorqual_ = {x: "True" for x in poor_qual}
    hdoublet_ = {x: "True" for x in h_doublet}
    ldoublet_ = {x: "True" for x in l_doublet}

    poorqual.update(poorqual_)
    hdoublet.update(hdoublet_)
    ldoublet.update(ldoublet_)

    adata_.obs["filter_contig_quality"] = pd.Series(poorqual)
    adata_.obs["filter_contig_VDJ"] = pd.Series(hdoublet)
    adata_.obs["filter_contig_VJ"] = pd.Series(ldoublet)

    drop_contig = list(set(flatten(drop_contig)))

    filter_ids = []
    if filter_contig:
        print("Finishing up filtering")
        if not filter_vj_chains:
            if filter_poorqualitycontig:
                filter_ids = list(set(h_doublet + poor_qual))
            else:
                filter_ids = list(set(h_doublet))
        else:
            if filter_poorqualitycontig:
                filter_ids = list(set(h_doublet + l_doublet + poor_qual))
            else:
                filter_ids = list(set(h_doublet + l_doublet))

        filter_ids = filter_ids + list(
            adata_[adata_.obs["filter_rna"].isin(TRUES)].obs_names
        )
        filter_ids = list(set(filter_ids))

        if filter_missing:
            dat = dat[dat["cell_id"].isin(adata_.obs_names)].copy()

        _dat = dat[~(dat["cell_id"].isin(filter_ids))].copy()
        _dat = _dat[~(_dat["sequence_id"].isin(drop_contig))].copy()

        # final check
        barcodes_final = list(set(_dat["cell_id"]))
        check_dat_barcodes = list(
            set(_dat[_dat["locus"].isin(HEAVYLONG)]["cell_id"])
        )
        filter_ids2 = list(set(barcodes_final) - set(check_dat_barcodes))
        _dat = _dat[~(_dat["cell_id"].isin(filter_ids2))].copy()

        if _dat.shape[0] == 0:
            raise IndexError(
                "No contigs passed filtering. Are you sure that the cell barcodes are matching?"
            )

        if os.path.isfile(str(data)):
            write_airr(
                _dat,
                "{}/{}_filtered.tsv".format(
                    os.path.dirname(data),
                    os.path.basename(data).split(".tsv")[0],
                ),
            )
        else:
            if save is not None:
                if save.endswith(".tsv"):
                    write_airr(_dat, str(save))
                else:
                    raise FileNotFoundError(
                        "{} not suitable. Please provide a file name that ends with .tsv".format(
                            str(save)
                        )
                    )
    else:
        _dat = dat.copy()

    if filter_contig:
        barcode1 = list(set(dat["cell_id"]))

    barcode2 = list(set(_dat["cell_id"]))

    if filter_contig:
        failed = list(set(barcode1) ^ set(barcode2))

    print("Initializing Dandelion object")
    out_dat = Dandelion(data=_dat, **kwargs)
    if data.__class__ == Dandelion:
        out_dat.germline = data.germline

    if adata_provided:
        bc_2 = {b: "True" for b in barcode2}
        if filter_contig:
            failed2 = {b: "False" for b in failed}
            bc_2.update(failed2)
        contig_check["contig_QC_pass"] = pd.Series(bc_2)
        contig_check.replace(np.nan, "No_contig", inplace=True)
        adata_.obs["contig_QC_pass"] = pd.Series(contig_check["contig_QC_pass"])
        adata_.obs["filter_contig"] = adata_.obs_names.isin(filter_ids)
        if filter_rna:
            # not saving the scanpy object because there's no need to at the moment
            out_adata = adata_[adata_.obs["filter_contig"].isin(FALSES)].copy()
        else:
            out_adata = adata_.copy()
        logg.info(
            " finished",
            time=start,
            deep=("Returning Dandelion and AnnData objects: \n"),
        )
        return (out_dat, out_adata)
    else:
        return out_dat


def quantify_mutations(
    self: Union[Dandelion, str, PathLike],
    split_locus: bool = False,
    sequence_column: Optional[str] = None,
    germline_column: Optional[str] = None,
    region_definition: Optional[str] = None,
    mutation_definition: Optional[str] = None,
    frequency: bool = False,
    combine: bool = True,
    **kwargs,
) -> Union[pd.DataFrame, Dandelion]:
    """
    Run basic mutation load analysis.

    Implemented in `shazam <https://shazam.readthedocs.io/en/stable/vignettes/Mutation-Vignette/>`__.

    Parameters
    ----------
    self : Dandelion, str, PathLike
        `Dandelion` object, file path to AIRR file.
    split_locus : bool
        whether to return the results for heavy chain and light chain separately. Default is False.
    sequence_column: str, Optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/observedMutations
    germline_column: str, Optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/observedMutations
    region_definition : str, Optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/IMGT_SCHEMES/
    mutation_definition : str, Optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/MUTATION_SCHEMES/
    frequency
        whether to return the results a frequency or counts. Default is True (frequency).
    combine
        whether to return the results for replacement and silent mutations separately (False). Default is True (sum).
    **kwargs
        passed to shazam::observedMutations.

    Returns
    -------
    `Dandelion` object with updated `.metadata` slot.
    """
    start = logg.info("Quantifying mutations")
    try:
        import rpy2
        from rpy2.robjects.packages import importr
        from rpy2.rinterface import NULL
        from rpy2.robjects import pandas2ri
    except:
        raise (
            ImportError(
                "Unable to initialise R instance. Please run this separately through R with Shazam's tutorial."
            )
        )

    sh = importr("shazam")
    base = importr("base")
    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    else:
        dat = load_data(self)

    pandas2ri.activate()
    warnings.filterwarnings("ignore")

    dat = sanitize_data(dat)

    if sequence_column is None:
        seq_ = "sequence_alignment"
    else:
        seq_ = sequence_column

    if germline_column is None:
        germline_ = "germline_alignment_d_mask"
    else:
        germline_ = germline_column

    if region_definition is None:
        reg_d = NULL
    else:

        reg_d = base.get(region_definition)

    if mutation_definition is None:
        mut_d = NULL
    else:
        mut_d = base.get(mutation_definition)

    if split_locus is False:
        dat = dat.where(dat.isna(), dat.astype(str))
        try:
            dat_r = pandas2ri.py2rpy(dat)
        except:
            dat = dat.astype(str)
            dat_r = pandas2ri.py2rpy(dat)

        results = sh.observedMutations(
            dat_r,
            sequenceColumn=seq_,
            germlineColumn=germline_,
            regionDefinition=reg_d,
            mutationDefinition=mut_d,
            frequency=frequency,
            combine=combine,
            **kwargs,
        )
        if rpy2.__version__ >= "3.4.5":
            from rpy2.robjects.conversion import localconverter

            with localconverter(
                rpy2.robjects.default_converter + pandas2ri.converter
            ):
                pd_df = rpy2.robjects.conversion.rpy2py(results)
        else:
            # pd_df = pandas2ri.rpy2py_dataframe(results)
            pd_df = results.copy()
    else:
        dat_h = dat[dat["locus"] == "IGH"]
        dat_l = dat[dat["locus"].isin(["IGK", "IGL"])]

        dat_h = dat_h.where(dat_h.isna(), dat_h.astype(str))
        try:
            dat_h_r = pandas2ri.py2rpy(dat_h)
        except:
            dat_h = dat_h.astype(str)
            dat_h_r = pandas2ri.py2rpy(dat_h)

        dat_l = dat_l.where(dat_l.isna(), dat_l.astype(str))
        try:
            dat_l_r = pandas2ri.py2rpy(dat_l)
        except:
            dat_l = dat_l.astype(str)
            dat_l_r = pandas2ri.py2rpy(dat_l)

        results_h = sh.observedMutations(
            dat_h_r,
            sequenceColumn=seq_,
            germlineColumn=germline_,
            regionDefinition=reg_d,
            mutationDefinition=mut_d,
            frequency=frequency,
            combine=combine,
            **kwargs,
        )
        results_l = sh.observedMutations(
            dat_l_r,
            sequenceColumn=seq_,
            germlineColumn=germline_,
            regionDefinition=reg_d,
            mutationDefinition=mut_d,
            frequency=frequency,
            combine=combine,
            **kwargs,
        )
        if rpy2.__version__ >= "3.4.5":
            from rpy2.robjects.conversion import localconverter

            with localconverter(
                rpy2.robjects.default_converter + pandas2ri.converter
            ):
                results_h = rpy2.robjects.conversion.rpy2py(results_h)
                results_l = rpy2.robjects.conversion.rpy2py(results_l)
        pd_df = pd.concat([results_h, results_l])

    pd_df.set_index("sequence_id", inplace=True, drop=False)
    # this doesn't actually catch overwritten columns
    cols_to_return = pd_df.columns.difference(dat.columns)
    if len(cols_to_return) < 1:
        cols_to_return = list(
            filter(re.compile("mu_.*").match, [c for c in pd_df.columns])
        )
    else:
        cols_to_return = cols_to_return

    res = {}
    if self.__class__ == Dandelion:
        for x in cols_to_return:
            res[x] = list(pd_df[x])
            # TODO: str will make it work for the back and forth conversion with rpy2. but maybe can use a better option
            self.data[x] = [str(r) for r in res[x]]
        self.data = sanitize_data(self.data)
        if split_locus is False:
            metadata_ = self.data[["cell_id"] + list(cols_to_return)]
        else:
            metadata_ = self.data[["locus", "cell_id"] + list(cols_to_return)]

        for x in cols_to_return:
            metadata_[x] = metadata_[x].astype(float)

        if split_locus is False:
            metadata_ = metadata_.groupby("cell_id").sum()
        else:
            metadata_ = metadata_.groupby(["locus", "cell_id"]).sum()
            metadatas = []
            for x in list(set(self.data["locus"])):
                tmp = metadata_.iloc[
                    metadata_.index.isin([x], level="locus"), :
                ]
                tmp.index = tmp.index.droplevel()
                tmp.columns = [c + "_" + str(x) for c in tmp.columns]
                metadatas.append(tmp)
            metadata_ = functools.reduce(
                lambda x, y: pd.merge(
                    x, y, left_index=True, right_index=True, how="outer"
                ),
                metadatas,
            )

        metadata_.index.name = None

        if self.metadata is None:
            self.metadata = metadata_
        else:
            for x in metadata_.columns:
                self.metadata[x] = pd.Series(metadata_[x])
        logg.info(
            " finished",
            time=start,
            deep=(
                "Updated Dandelion object: \n"
                "   'data', contig-indexed clone table\n"
                "   'metadata', cell-indexed clone table\n"
            ),
        )
    else:
        for x in cols_to_return:
            res[x] = list(pd_df[x])
            # TODO: str will make it work for the back and forth conversion with rpy2. but maybe can use a better option
            dat[x] = [str(r) for r in res[x]]
        # dat = sanitize_data(dat)
        if self.__class__ == pd.DataFrame:
            logg.info(" finished", time=start, deep=("Returning DataFrame\n"))
            return dat
        elif os.path.isfile(self):
            logg.info(
                " finished",
                time=start,
                deep=("saving DataFrame at {}\n".format(str(self))),
            )
            write_airr(dat, self)


def calculate_threshold(
    self: Union[Dandelion, pd.DataFrame, str],
    mode: Literal["single-cell", "heavy"] = "single-cell",
    manual_threshold: Optional[float] = None,
    VJthenLen: bool = False,
    onlyHeavy: bool = False,
    model: Optional[
        Literal[
            "ham",
            "aa",
            "hh_s1f",
            "hh_s5f",
            "mk_rs1nf",
            "hs1f_compat",
            "m1n_compat",
        ]
    ] = None,
    normalize_method: Optional[Literal["len"]] = None,
    threshold_method: Optional[Literal["gmm", "density"]] = None,
    edge: Optional[float] = None,
    cross: Optional[Sequence] = None,
    subsample: Optional[int] = None,
    threshold_model: Optional[
        Literal["norm-norm", "norm-gamma", "gamma-norm", "gamma-gamma"]
    ] = None,
    cutoff: Optional[Literal["optimal", "intersect", "user"]] = None,
    sensitivity: Optional[float] = None,
    specificity: Optional[float] = None,
    plot: bool = True,
    plot_group: Optional[str] = None,
    figsize: Tuple[Union[int, float], Union[int, float]] = (4.5, 2.5),
    ncpu: int = 1,
    **kwargs,
) -> Dandelion:
    """
    Calculating nearest neighbor distances for tuning clonal assignment with `shazam`.

    <https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/>`__.
    Runs the following:

    distToNearest
        Get non-zero distance of every heavy chain (IGH) sequence (as defined by sequenceColumn) to its nearest sequence
        in a partition of heavy chains sharing the same V gene, J gene, and junction length (VJL), or in a partition of
        single cells with heavy chains sharing the same heavy chain VJL combination, or of single cells with heavy and
        light chains sharing the same heavy chain VJL and light chain VJL combinations.
    findThreshold
        automtically determines an optimal threshold for clonal assignment of Ig sequences using a vector of nearest
        neighbor distances. It provides two alternative methods using either a Gamma/Gaussian Mixture Model fit
        (threshold_method="gmm") or kernel density fit (threshold_method="density").

    Parameters
    ----------
    self : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones
        have been determined.
    mode : Literal, str
        accepts one of "heavy" or "single-cell".
        Refer to https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette.
    manual_threshold : float, Optional
        value to manually plot in histogram.
    VJthenLen : bool
        logical value specifying whether to perform partitioning as a 2-stage process.
        If True, partitions are made first based on V and J gene, and then further split
        based on junction lengths corresponding to sequenceColumn.
        If False, perform partition as a 1-stage process during which V gene, J gene, and junction length
        are used to create partitions simultaneously.
        Defaults to False.
    onlyHeavy : bool
        use only the IGH (BCR) or TRB/TRD (TCR) sequences for grouping. Only applicable to single-cell mode.
        See groupGenes for further details.
    model : str, Optional
        underlying SHM model, which must be one of "ham","aa","hh_s1f","hh_s5f","mk_rs1nf","hs1f_compat","m1n_compat".
    normalize_method : str, Optional
        method of normalization. The default is "len", which divides the distance by the length of the sequence group.
        If "none" then no normalization if performed.
    threshold_method : str, Optional
        string defining the method to use for determining the optimal threshold. One of "gmm" or "density".
    edge : float, Optional
        upper range as a fraction of the data density to rule initialization of Gaussian fit parameters.
        Default value is 0.9 (or 90). Applies only when threshold_method="density".
    cross : Sequence, Optional
        supplementary nearest neighbor distance vector output from distToNearest for initialization of the Gaussian fit
        parameters. Applies only when method="gmm".
    subsample : int, Optional
        maximum number of distances to subsample to before threshold detection.
    threshold_model : str, Optional
        allows the user to choose among four possible combinations of fitting curves: "norm-norm", "norm-gamma",
        "gamma-norm", and "gamma-gamma". Applies only when method="gmm".
    cutoff : str, Optional
        method to use for threshold selection: the optimal threshold "optimal", the intersection point of the two fitted
        curves "intersect", or a value defined by user for one of the sensitivity or specificity "user". Applies only
        when method="gmm".
    sensitivity : float, Optional
        sensitivity required. Applies only when method="gmm" and cutoff="user".
    specificity : float, Optional
        specificity required. Applies only when method="gmm" and cutoff="user".
    plot : bool
        whether or not to return plot.
    plot_group : str, Optional
        determines the fill color and facets.
    figsize : Tuple[Union[int,float], Union[int,float]]
        size of plot. Default is (4.5, 2.5).
    ncpu : float
        number of cpus to run `distToNearest`. defaults to 1.
    **kwargs
        passed to shazam's `distToNearest <https://shazam.readthedocs.io/en/stable/topics/distToNearest/>`__.

    Returns
    -------
        `Dandelion` object object with distance threshold value in `.threshold`.

        If plot = True,plotnine plot showing histogram of length normalized ham model distance threshold.
    """
    start = logg.info("Calculating threshold")
    try:
        import rpy2
        from rpy2.robjects.packages import importr
        from rpy2.rinterface import NULL
        from rpy2.robjects import pandas2ri, FloatVector
    except:
        raise (
            ImportError(
                "Unable to initialise R instance. Please run this separately through R with Shazam's tutorial."
            )
        )

    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    elif self.__class__ == pd.DataFrame or os.path.isfile(str(self)):
        dat = load_data(self)
        warnings.filterwarnings("ignore")

    sh = importr("shazam")
    pandas2ri.activate()
    if "v_call_genotyped" in dat.columns:
        v_call = "v_call_genotyped"
    else:
        v_call = "v_call"
    if model is None:
        model_ = "ham"
    else:
        model_ = model
    if normalize_method is None:
        norm_ = "len"
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

    if mode == "heavy":
        dat_h = dat[dat["locus"].isin(["IGH", "TRB", "TRD"])].copy()
        try:
            dat_h_r = pandas2ri.py2rpy(dat_h)
        except:
            dat_h = dat_h.astype(str)
            dat_h_r = pandas2ri.py2rpy(dat_h)

        dist_ham = sh.distToNearest(
            dat_h_r, vCallColumn=v_call, model=model_, normalize=norm_, **kwargs
        )
    elif mode == "single-cell":
        try:
            dat_r = pandas2ri.py2rpy(dat)
        except:
            dat = dat.astype(str)
            dat_r = pandas2ri.py2rpy(dat)
        try:
            dist_ham = sh.distToNearest(
                dat_r,
                cellIdColumn="cell_id",
                locusColumn="locus",
                VJthenLen=VJthenLen,
                vCallColumn=v_call,
                onlyHeavy=onlyHeavy,
                normalize=norm_,
                model=model_,
                nproc=ncpu,
                **kwargs,
            )
        except:
            print(
                "Rerun this after filtering. For now, switching to heavy mode."
            )
            dat_h = dat[dat["locus"].isin(["IGH", "TRB", "TRD"])].copy()
            try:
                dat_h_r = pandas2ri.py2rpy(dat_h)
            except:
                dat_h = dat_h.astype(str)
                dat_h_r = pandas2ri.py2rpy(dat_h)

            dist_ham = sh.distToNearest(
                dat_h_r,
                vCallColumn=v_call,
                model=model_,
                normalize=norm_,
                nproc=ncpu,
                **kwargs,
            )
    if rpy2.__version__ >= "3.4.5":
        from rpy2.robjects.conversion import localconverter

        with localconverter(
            rpy2.robjects.default_converter + pandas2ri.converter
        ):
            dist_ham = rpy2.robjects.conversion.rpy2py(dist_ham)
    # Find threshold using density method
    dist = np.array(dist_ham["dist_nearest"])
    if threshold_method_ == "density":
        if edge is None:
            edge_ = 0.9
        else:
            edge_ = edge
        dist_threshold = sh.findThreshold(
            FloatVector(dist[~np.isnan(dist)]),
            method=threshold_method_,
            subsample=subsample_,
            edge=edge_,
        )
        threshold = np.array(dist_threshold.slots["threshold"])[0]
        if np.isnan(threshold):
            print(
                "      Threshold method 'density' did not return with any values. Switching to method = 'gmm'."
            )
            threshold_method_ = "gmm"
            if threshold_model is None:
                threshold_model_ = "gamma-gamma"
            else:
                threshold_model_ = threshold_model
            if cross is None:
                cross_ = NULL
            else:
                cross_ = cross
            if cutoff is None:
                cutoff_ = "optimal"
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
            dist_threshold = sh.findThreshold(
                FloatVector(dist[~np.isnan(dist)]),
                method=threshold_method_,
                model=threshold_model_,
                cross=cross_,
                subsample=subsample_,
                cutoff=cutoff_,
                sen=sen_,
                spc=spc_,
            )
            if rpy2.__version__ >= "3.4.5":
                from rpy2.robjects.conversion import localconverter

                with localconverter(
                    rpy2.robjects.default_converter + pandas2ri.converter
                ):
                    dist_threshold = rpy2.robjects.conversion.rpy2py(
                        dist_threshold
                    )

            threshold = np.array(dist_threshold.slots["threshold"])[0]
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
            cutoff_ = "optimal"
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
        dist_threshold = sh.findThreshold(
            FloatVector(dist[~np.isnan(dist)]),
            method=threshold_method_,
            model=threshold_model_,
            cross=cross_,
            subsample=subsample_,
            cutoff=cutoff_,
            sen=sen_,
            spc=spc_,
        )
        if rpy2.__version__ >= "3.4.5":
            from rpy2.robjects.conversion import localconverter

            with localconverter(
                rpy2.robjects.default_converter + pandas2ri.converter
            ):
                dist_threshold = rpy2.robjects.conversion.rpy2py(dist_threshold)
        threshold = np.array(dist_threshold.slots["threshold"])[0]
    if np.isnan(threshold):
        raise ValueError(
            "Automatic thresholding failed. Please visually inspect the resulting distribution fits"
            + " and choose a threshold value manually."
        )
    # dist_ham = pandas2ri.rpy2py_dataframe(dist_ham)

    if manual_threshold is None:
        tr = threshold
    else:
        tr = manual_threshold

    if plot:
        options.figure_size = figsize
        if plot_group is None:
            plot_group = "sample_id"
        else:
            plot_group = plot_group

        print(
            (
                ggplot(dist_ham, aes("dist_nearest", fill=str(plot_group)))
                + theme_bw()
                + xlab("Grouped Hamming distance")
                + ylab("Count")
                + geom_histogram(binwidth=0.01)
                + geom_vline(
                    xintercept=tr, linetype="dashed", color="blue", size=0.5
                )
                + annotate(
                    "text",
                    x=tr + 0.02,
                    y=10,
                    label="Threshold:\n" + str(np.around(tr, decimals=2)),
                    size=8,
                    color="Blue",
                )
                + facet_wrap("~" + str(plot_group), scales="free_y")
                + theme(legend_position="none")
            )
        )
    else:
        print(
            "Automatic Threshold : "
            + str(np.around(threshold, decimals=2))
            + "\n method = "
            + str(threshold_method_)
        )
    if self.__class__ == Dandelion:
        self.threshold = tr
        logg.info(
            " finished",
            time=start,
            deep=(
                "Updated Dandelion object: \n"
                "   'threshold', threshold value for tuning clonal assignment\n"
            ),
        )
    else:
        output = Dandelion(dat)
        output.threshold = tr
        return output


class FilterContigs:
    """
    `FilterContigs` class object.

    Main class object to run filter_contigs.

    """

    def __init__(
        self,
        data,
        keep_highest_umi,
        umi_foldchange_cutoff,
        filter_poorqualitycontig,
    ):
        self.Cell = Tree()
        self.poor_qual = []
        self.h_doublet = []
        self.l_doublet = []
        self.drop_contig = []
        self.umi_adjustment = {}
        if "v_call_genotyped" in data.columns:
            v_dict = dict(zip(data["sequence_id"], data["v_call_genotyped"]))
        else:
            v_dict = dict(zip(data["sequence_id"], data["v_call"]))
        d_dict = dict(zip(data["sequence_id"], data["d_call"]))
        j_dict = dict(zip(data["sequence_id"], data["j_call"]))
        c_dict = dict(zip(data["sequence_id"], data["c_call"]))
        for contig, row in tqdm(data.iterrows(), desc="Preparing data"):
            cell = row["cell_id"]
            if row["locus"] in HEAVYLONG:
                if row["productive"] in TRUES:
                    self.Cell[cell]["VDJ"]["P"][contig].update(row)
                elif row["productive"] in FALSES:
                    self.Cell[cell]["VDJ"]["NP"][contig].update(row)
            elif row["locus"] in LIGHTSHORT:
                if row["productive"] in TRUES:
                    self.Cell[cell]["VJ"]["P"][contig].update(row)
                elif row["productive"] in FALSES:
                    self.Cell[cell]["VJ"]["NP"][contig].update(row)
        for cell in tqdm(
            self.Cell, desc="Scanning for poor quality/ambiguous contigs"
        ):
            if len(self.Cell[cell]["VDJ"]["P"]) > 0:
                data1 = pd.DataFrame(
                    [
                        self.Cell[cell]["VDJ"]["P"][x]
                        for x in self.Cell[cell]["VDJ"]["P"]
                        if isinstance(self.Cell[cell]["VDJ"]["P"][x], dict)
                    ],
                    index=[
                        self.Cell[cell]["VDJ"]["P"][x]["sequence_id"]
                        for x in self.Cell[cell]["VDJ"]["P"]
                        if isinstance(self.Cell[cell]["VDJ"]["P"][x], dict)
                    ],
                )
                h_p = list(data1["sequence_id"])
                h_umi_p = [
                    int(x) for x in pd.to_numeric(data1["duplicate_count"])
                ]
                h_ccall_p = list(data1["c_call"])
                if len(h_p) > 1:
                    if "sequence_alignment" in data1:
                        h_seq_p = list(data1["sequence_alignment"])
                        if len(set(h_seq_p)) == 1:
                            if len(set(h_ccall_p)) == 1:
                                highest_umi_h = max(h_umi_p)
                                highest_umi_h_idx = [
                                    i
                                    for i, j in enumerate(h_umi_p)
                                    if j == highest_umi_h
                                ]
                                keep_index_h = highest_umi_h_idx[0]
                                self.drop_contig.append(
                                    h_p[:keep_index_h] + h_p[keep_index_h:]
                                )
                                keep_hc_contig = h_p[keep_index_h]
                                data1[keep_hc_contig, "duplicate_count"] = int(
                                    np.sum(
                                        h_umi_p[:keep_index_h]
                                        + h_umi_p[keep_index_h:]
                                    )
                                )
                                self.umi_adjustment.update(
                                    {
                                        keep_hc_contig: int(
                                            np.sum(
                                                h_umi_p[:keep_index_h]
                                                + h_umi_p[keep_index_h:]
                                            )
                                        )
                                    }
                                )
                                # refresh
                                data1 = pd.DataFrame(
                                    [data1.loc[keep_hc_contig]]
                                )
                                h_p = list(data1["sequence_id"])
                                h_umi_p = [
                                    int(x)
                                    for x in pd.to_numeric(
                                        data1["duplicate_count"]
                                    )
                                ]
                                h_ccall_p = list(data1["c_call"])
                    if len(h_p) > 1:
                        highest_umi_h = max(h_umi_p)
                        highest_umi_idx = [
                            i
                            for i, j in enumerate(h_umi_p)
                            if j == highest_umi_h
                        ]
                        keep_index_h = highest_umi_idx[0]
                        keep_hc_contig = h_p[keep_index_h]
                        umi_test = [
                            int(highest_umi_h) / x < umi_foldchange_cutoff
                            for x in h_umi_p[:keep_index_h]
                            + h_umi_p[keep_index_h:]
                        ]
                        sum_umi = sum(h_umi_p)
                        if "IGHM" and "IGHD" in h_ccall_p:
                            if all(
                                cc_ == "IGHM" or cc_ == "IGHD"
                                for cc_ in h_ccall_p
                            ):
                                pass
                        else:
                            if len(highest_umi_idx) > 1:
                                self.h_doublet.append(cell)
                            if sum_umi < 4:
                                self.h_doublet.append(cell)
                            if any(umi_test):
                                self.h_doublet.append(cell)
                            if len(highest_umi_idx) == 1:
                                other_umi_idx = [
                                    i
                                    for i, j in enumerate(h_umi_p)
                                    if j != highest_umi_h
                                ]
                                umi_test_ = [
                                    highest_umi_h / x >= umi_foldchange_cutoff
                                    for x in h_umi_p[:keep_index_h]
                                    + h_umi_p[keep_index_h:]
                                ]
                                umi_test_dict = dict(
                                    zip(other_umi_idx, umi_test_)
                                )
                                for otherindex in umi_test_dict:
                                    if umi_test_dict[otherindex]:
                                        if keep_highest_umi:
                                            self.drop_contig.append(
                                                h_p[otherindex]
                                            )
                                            # refresh
                                data1 = pd.DataFrame(
                                    [data1.loc[keep_hc_contig]]
                                )
                                h_p = list(data1["sequence_id"])
            if len(self.Cell[cell]["VDJ"]["NP"]) > 0:
                data2 = pd.DataFrame(
                    [
                        self.Cell[cell]["VDJ"]["NP"][x]
                        for x in self.Cell[cell]["VDJ"]["NP"]
                        if isinstance(self.Cell[cell]["VDJ"]["NP"][x], dict)
                    ],
                    index=[
                        self.Cell[cell]["VDJ"]["NP"][x]["sequence_id"]
                        for x in self.Cell[cell]["VDJ"]["NP"]
                        if isinstance(self.Cell[cell]["VDJ"]["NP"][x], dict)
                    ],
                )
                h_np = list(data2["sequence_id"])
                h_umi_np = [
                    int(x) for x in pd.to_numeric(data2["duplicate_count"])
                ]
                if len(h_np) > 1:
                    highest_umi_h = max(h_umi_np)
                    highest_umi_idx = [
                        i for i, j in enumerate(h_umi_np) if j == highest_umi_h
                    ]
                    if len(highest_umi_idx) == 1:
                        keep_index_h = highest_umi_idx[0]
                        keep_hc_contig = h_np[keep_index_h]
                        other_umi_idx = [
                            i
                            for i, j in enumerate(h_umi_np)
                            if j != highest_umi_h
                        ]
                        umi_test_ = [
                            highest_umi_h / x >= umi_foldchange_cutoff
                            for x in h_umi_np[:keep_index_h]
                            + h_umi_np[keep_index_h:]
                        ]
                        umi_test_dict = dict(zip(other_umi_idx, umi_test_))
                        for otherindex in umi_test_dict:
                            if umi_test_dict[otherindex]:
                                self.drop_contig.append(h_np[otherindex])
                        # refresh
                        data2 = pd.DataFrame([data2.loc[keep_hc_contig]])
                        h_np = list(data2["sequence_id"])
                        h_umi_np = [
                            int(x)
                            for x in pd.to_numeric(data2["duplicate_count"])
                        ]
            if len(self.Cell[cell]["VJ"]["P"]) > 0:
                data3 = pd.DataFrame(
                    [
                        self.Cell[cell]["VJ"]["P"][x]
                        for x in self.Cell[cell]["VJ"]["P"]
                        if isinstance(self.Cell[cell]["VJ"]["P"][x], dict)
                    ],
                    index=[
                        self.Cell[cell]["VJ"]["P"][x]["sequence_id"]
                        for x in self.Cell[cell]["VJ"]["P"]
                        if isinstance(self.Cell[cell]["VJ"]["P"][x], dict)
                    ],
                )
                l_p = list(data3["sequence_id"])
                l_umi_p = [
                    int(x) for x in pd.to_numeric(data3["duplicate_count"])
                ]
                if len(l_p) > 1:
                    if "sequence_alignment" in data3:
                        l_seq_p = list(data3["sequence_alignment"])
                        if len(list(set(l_seq_p))) == 1:
                            highest_umi_l = max(l_umi_p)
                            highest_umi_l_idx = [
                                i
                                for i, j in enumerate(l_umi_p)
                                if j == highest_umi_l
                            ]
                            keep_index_l = highest_umi_l_idx[0]
                            self.drop_contig.append(
                                l_p[:keep_index_l] + l_p[keep_index_l:]
                            )
                            keep_lc_contig = l_p[keep_index_l]
                            data3.at[keep_lc_contig, "duplicate_count"] = int(
                                np.sum(
                                    l_umi_p[:keep_index_l]
                                    + l_umi_p[keep_index_l:]
                                )
                            )
                            self.umi_adjustment.update(
                                {
                                    keep_lc_contig: int(
                                        np.sum(
                                            l_umi_p[:keep_index_l]
                                            + l_umi_p[keep_index_l:]
                                        )
                                    )
                                }
                            )
                            # refresh
                            data3 = pd.DataFrame([data3.loc[keep_lc_contig]])
                            l_p = list(data3["sequence_id"])
                            l_umi_p = [
                                int(x)
                                for x in pd.to_numeric(data3["duplicate_count"])
                            ]
                    if len(l_p) > 1:
                        highest_umi_l = max(l_umi_p)
                        highest_umi_l_idx = [
                            i
                            for i, j in enumerate(l_umi_p)
                            if j == highest_umi_l
                        ]
                        keep_index_l = highest_umi_l_idx[0]
                        keep_lc_contig = l_p[keep_index_l]
                        umi_test = [
                            highest_umi_l / x < umi_foldchange_cutoff
                            for x in l_umi_p[:keep_index_l]
                            + l_umi_p[keep_index_l:]
                        ]
                        sum_umi = sum(l_umi_p)
                        if len(highest_umi_l_idx) > 1:
                            self.l_doublet.append(cell)
                        if sum_umi < 4:
                            self.l_doublet.append(cell)
                        if any(umi_test):
                            self.l_doublet.append(cell)
                        if len(highest_umi_l_idx) == 1:
                            other_umi_idx_l = [
                                i
                                for i, j in enumerate(l_umi_p)
                                if j != highest_umi_l
                            ]
                            umi_test_l = [
                                highest_umi_l / x >= umi_foldchange_cutoff
                                for x in l_umi_p[:keep_index_l]
                                + l_umi_p[keep_index_l:]
                            ]
                            umi_test_dict_l = dict(
                                zip(other_umi_idx_l, umi_test_l)
                            )
                            for otherindex in umi_test_dict_l:
                                if umi_test_dict_l[otherindex]:
                                    if keep_highest_umi:
                                        self.drop_contig.append(l_p[otherindex])
                                        # refresh
                            data3 = pd.DataFrame([data3.loc[keep_lc_contig]])
                            l_p = list(data3["sequence_id"])
            if len(self.Cell[cell]["VJ"]["NP"]) > 0:
                data4 = pd.DataFrame(
                    [
                        self.Cell[cell]["VJ"]["NP"][x]
                        for x in self.Cell[cell]["VJ"]["NP"]
                        if isinstance(self.Cell[cell]["VJ"]["NP"][x], dict)
                    ],
                    index=[
                        self.Cell[cell]["VJ"]["NP"][x]["sequence_id"]
                        for x in self.Cell[cell]["VJ"]["NP"]
                        if isinstance(self.Cell[cell]["VJ"]["NP"][x], dict)
                    ],
                )
                l_np = list(data4["sequence_id"])
                l_umi_np = [
                    int(x) for x in pd.to_numeric(data4["duplicate_count"])
                ]
                if len(l_np) > 1:
                    highest_umi_l = max(l_umi_np)
                    highest_umi_l_idx = [
                        i for i, j in enumerate(l_umi_np) if j == highest_umi_l
                    ]
                    keep_index_l = highest_umi_l_idx[0]
                    keep_lc_contig = l_np[keep_index_l]
                    other_umi_idx_l = [
                        i for i, j in enumerate(l_umi_np) if j != highest_umi_l
                    ]
                    umi_test_l = [
                        highest_umi_l / x >= umi_foldchange_cutoff
                        for x in l_umi_np[:keep_index_l]
                        + l_umi_np[keep_index_l:]
                    ]
                    if len(highest_umi_l_idx) == 1:
                        umi_test_dict_l = dict(zip(other_umi_idx_l, umi_test_l))
                        for otherindex in umi_test_dict_l:
                            if umi_test_dict_l[otherindex]:
                                if keep_highest_umi:
                                    self.drop_contig.append(l_np[otherindex])
                        data4 = pd.DataFrame([data4.loc[keep_lc_contig]])
                        l_np = list(data4["sequence_id"])

            if "h_p" not in locals():
                h_p = []
            if "l_p" not in locals():
                l_p = []
            if "h_np" not in locals():
                h_np = []
            if "l_np" not in locals():
                l_np = []

            # marking doublets defined by VJ chains
            if (len(h_p) == 1) & (len(l_p) > 1):
                self.l_doublet.append(cell)

            # marking poor bcr quality, defined as those with only VJ chains, those
            # that were have conflicting assignment of locus and V(D)J v-, d-, j- and c- calls,
            # and also those that are missing j calls (to catch non-productive).
            if len(h_p) < 1:
                if filter_poorqualitycontig:
                    self.poor_qual.append(cell)
                self.drop_contig.append(l_p)
            if len(h_p) == 1:
                v = v_dict[h_p[0]]
                j = j_dict[h_p[0]]
                d = d_dict[h_p[0]]
                c = c_dict[h_p[0]]
                if present(v):
                    if not re.search("IGH|TR[BD]|TRAV.*/DV", v):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(cell)
                        self.drop_contig.append(l_p)
                        self.drop_contig.append(h_p)
                if present(d):
                    if not re.search("IGH|TR[BD]", d):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(cell)
                        self.drop_contig.append(l_p)
                        self.drop_contig.append(h_p)
                if present(j):
                    if not re.search("IGH|TR[BD]", j):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(cell)
                        self.drop_contig.append(l_p)
                        self.drop_contig.append(h_p)
                if present(c):
                    if not re.search("IGH|TR[BD]", c):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(cell)
                        self.drop_contig.append(l_p)
                        self.drop_contig.append(h_p)

                if present(j):
                    if present(v):
                        if not_same_call(v, j, "IGH"):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(l_p)
                            self.drop_contig.append(h_p)
                        elif not_same_call(v, j, "TRB"):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(l_p)
                            self.drop_contig.append(h_p)
                        elif not_same_call(v, j, "TRD"):
                            if not re.search("TRAV.*/DV", v):
                                if filter_poorqualitycontig:
                                    self.poor_qual.append(cell)
                                self.drop_contig.append(l_p)
                                self.drop_contig.append(h_p)

                    if present(d):
                        if not_same_call(d, j, "IGH"):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(l_p)
                            self.drop_contig.append(h_p)
                        elif not_same_call(d, j, "TRB"):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(l_p)
                            self.drop_contig.append(h_p)
                        elif not_same_call(d, j, "TRD"):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(l_p)
                            self.drop_contig.append(h_p)
                else:
                    if filter_poorqualitycontig:
                        self.poor_qual.append(cell)
                    self.drop_contig.append(l_p)
                    self.drop_contig.append(h_p)

            if len(h_p) > 1:
                for hx in h_p:
                    v = v_dict[hx]
                    d = d_dict[hx]
                    j = j_dict[hx]
                    c = c_dict[hx]
                    if present(v):
                        if not re.search("IGH|TR[BD]|TRAV.*/DV", v):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(hx)
                    if present(d):
                        if not re.search("IGH|TR[BD]", d):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(hx)
                    if present(j):
                        if not re.search("IGH|TR[BD]", j):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(hx)
                    if present(c):
                        if not re.search("IGH|TR[BD]", c):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(hx)
                    if present(j):
                        if present(v):
                            if not_same_call(v, j, "IGH"):
                                if filter_poorqualitycontig:
                                    self.poor_qual.append(cell)
                                self.drop_contig.append(hx)
                            elif not_same_call(v, j, "TRB"):
                                if filter_poorqualitycontig:
                                    self.poor_qual.append(cell)
                                self.drop_contig.append(hx)
                            elif not_same_call(v, j, "TRD"):
                                if not re.search("TRAV.*/DV", v):
                                    if filter_poorqualitycontig:
                                        self.poor_qual.append(cell)
                                    self.drop_contig.append(hx)
                        if present(d):
                            if not_same_call(d, j, "IGH"):
                                if filter_poorqualitycontig:
                                    self.poor_qual.append(cell)
                                self.drop_contig.append(hx)
                            elif not_same_call(d, j, "TRB"):
                                if filter_poorqualitycontig:
                                    self.poor_qual.append(cell)
                                self.drop_contig.append(hx)
                            elif not_same_call(d, j, "TRD"):
                                if filter_poorqualitycontig:
                                    self.poor_qual.append(cell)
                                self.drop_contig.append(hx)
                    else:
                        if filter_poorqualitycontig:
                            self.poor_qual.append(cell)
                        self.drop_contig.append(hx)

            if len(h_np) > 0:
                for hx in h_np:
                    v = v_dict[hx]
                    d = d_dict[hx]
                    j = j_dict[hx]
                    c = c_dict[hx]
                    if present(v):
                        if not re.search("IGH|TR[BD]|TRAV.*/DV", v):
                            self.drop_contig.append(hx)
                    if present(d):
                        if not re.search("IGH|TR[BD]", d):
                            self.drop_contig.append(hx)
                    if present(j):
                        if not re.search("IGH|TR[BD]", j):
                            self.drop_contig.append(hx)
                    if present(c):
                        if not re.search("IGH|TR[BD]", c):
                            self.drop_contig.append(hx)

                    if present(j):
                        if present(v):
                            if not_same_call(v, j, "IGH"):
                                self.drop_contig.append(hx)
                            elif not_same_call(v, j, "TRB"):
                                self.drop_contig.append(hx)
                            elif not_same_call(v, j, "TRD"):
                                if not re.search("TRAV.*/DV", v):
                                    self.drop_contig.append(hx)
                        if present(d):
                            if not_same_call(d, j, "IGH"):
                                self.drop_contig.append(hx)
                            elif not_same_call(d, j, "TRB"):
                                self.drop_contig.append(hx)
                            elif not_same_call(d, j, "TRD"):
                                self.drop_contig.append(hx)
                    else:
                        self.drop_contig.append(hx)
            if len(l_p) > 0:
                for lx in l_p:
                    v = v_dict[lx]
                    j = j_dict[lx]
                    c = c_dict[lx]
                    if present(v):
                        if re.search("IGH|TR[BD]", v):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(lx)
                    if present(j):
                        if re.search("IGH|TR[BD]", j):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(lx)
                    if present(c):
                        if re.search("IGH|TR[BD]", c):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(cell)
                            self.drop_contig.append(lx)

                    if present(j):
                        if present(v):
                            if not_same_call(v, j, "IGK"):
                                if filter_poorqualitycontig:
                                    self.poor_qual.append(cell)
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "IGL"):
                                if filter_poorqualitycontig:
                                    self.poor_qual.append(cell)
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "TRA"):
                                if filter_poorqualitycontig:
                                    self.poor_qual.append(cell)
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "TRG"):
                                if filter_poorqualitycontig:
                                    self.poor_qual.append(cell)
                                self.drop_contig.append(lx)
                    else:
                        if filter_poorqualitycontig:
                            self.poor_qual.append(cell)
                        self.drop_contig.append(lx)

            if len(l_np) > 0:
                for lx in l_np:
                    v = v_dict[lx]
                    j = j_dict[lx]
                    c = c_dict[lx]
                    if present(v):
                        if re.search("IGH|TR[BD]", v):
                            self.drop_contig.append(lx)
                    if present(j):
                        if re.search("IGH|TR[BD]", j):
                            self.drop_contig.append(lx)
                    if present(c):
                        if re.search("IGH|TR[BD]", c):
                            self.drop_contig.append(lx)

                    if present(j):
                        if present(v):
                            if not_same_call(v, j, "IGK"):
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "IGL"):
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "TRA"):
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "TRG"):
                                self.drop_contig.append(lx)
                    else:
                        self.drop_contig.append(lx)


class FilterContigsLite:
    """
    `FilterContigsLite` class object.

    Main class object to run filter_contigs, lite mode.

    """

    def __init__(self, data):
        self.Cell = Tree()
        self.poor_qual = []
        self.h_doublet = []
        self.l_doublet = []
        self.drop_contig = []
        self.umi_adjustment = {}
        if "v_call_genotyped" in data.columns:
            v_dict = dict(zip(data["sequence_id"], data["v_call_genotyped"]))
        else:
            v_dict = dict(zip(data["sequence_id"], data["v_call"]))
        d_dict = dict(zip(data["sequence_id"], data["d_call"]))
        j_dict = dict(zip(data["sequence_id"], data["j_call"]))
        c_dict = dict(zip(data["sequence_id"], data["c_call"]))
        for contig, row in tqdm(data.iterrows(), desc="Preparing data"):
            cell = row["cell_id"]
            if row["locus"] in HEAVYLONG:
                if row["productive"] in TRUES:
                    self.Cell[cell]["VDJ"]["P"][contig].update(row)
                elif row["productive"] in FALSES:
                    self.Cell[cell]["VDJ"]["NP"][contig].update(row)
            elif row["locus"] in LIGHTSHORT:
                if row["productive"] in TRUES:
                    self.Cell[cell]["VJ"]["P"][contig].update(row)
                elif row["productive"] in FALSES:
                    self.Cell[cell]["VJ"]["NP"][contig].update(row)
        for cell in tqdm(
            self.Cell, desc="Scanning for poor quality/ambiguous contigs"
        ):
            if len(self.Cell[cell]["VDJ"]["P"]) > 0:
                data1 = pd.DataFrame(
                    [
                        self.Cell[cell]["VDJ"]["P"][x]
                        for x in self.Cell[cell]["VDJ"]["P"]
                        if isinstance(self.Cell[cell]["VDJ"]["P"][x], dict)
                    ],
                    index=[
                        self.Cell[cell]["VDJ"]["P"][x]["sequence_id"]
                        for x in self.Cell[cell]["VDJ"]["P"]
                        if isinstance(self.Cell[cell]["VDJ"]["P"][x], dict)
                    ],
                )
                h_p = list(data1["sequence_id"])
                h_umi_p = [
                    int(x) for x in pd.to_numeric(data1["duplicate_count"])
                ]
                h_ccall_p = list(data1["c_call"])
                if len(h_p) > 1:
                    if "sequence_alignment" in data1:
                        h_seq_p = list(data1["sequence_alignment"])
                        if len(set(h_seq_p)) == 1:
                            if len(set(h_ccall_p)) == 1:
                                highest_umi_h = max(h_umi_p)
                                highest_umi_h_idx = [
                                    i
                                    for i, j in enumerate(h_umi_p)
                                    if j == highest_umi_h
                                ]
                                keep_index_h = highest_umi_h_idx[0]
                                self.drop_contig.append(
                                    h_p[:keep_index_h] + h_p[keep_index_h:]
                                )
                                keep_hc_contig = h_p[keep_index_h]
                                data1[keep_hc_contig, "duplicate_count"] = int(
                                    np.sum(
                                        h_umi_p[:keep_index_h]
                                        + h_umi_p[keep_index_h:]
                                    )
                                )
                                self.umi_adjustment.update(
                                    {
                                        keep_hc_contig: int(
                                            np.sum(
                                                h_umi_p[:keep_index_h]
                                                + h_umi_p[keep_index_h:]
                                            )
                                        )
                                    }
                                )
                                # refresh
                                data1 = pd.DataFrame(
                                    [data1.loc[keep_hc_contig]]
                                )
                                h_p = list(data1["sequence_id"])
                                h_umi_p = [
                                    int(x)
                                    for x in pd.to_numeric(
                                        data1["duplicate_count"]
                                    )
                                ]
            if len(self.Cell[cell]["VDJ"]["NP"]) > 0:
                data2 = pd.DataFrame(
                    [
                        self.Cell[cell]["VDJ"]["NP"][x]
                        for x in self.Cell[cell]["VDJ"]["NP"]
                        if isinstance(self.Cell[cell]["VDJ"]["NP"][x], dict)
                    ],
                    index=[
                        self.Cell[cell]["VDJ"]["NP"][x]["sequence_id"]
                        for x in self.Cell[cell]["VDJ"]["NP"]
                        if isinstance(self.Cell[cell]["VDJ"]["NP"][x], dict)
                    ],
                )
                h_np = list(data2["sequence_id"])
                h_umi_np = [
                    int(x) for x in pd.to_numeric(data2["duplicate_count"])
                ]
            if len(self.Cell[cell]["VJ"]["P"]) > 0:
                data3 = pd.DataFrame(
                    [
                        self.Cell[cell]["VJ"]["P"][x]
                        for x in self.Cell[cell]["VJ"]["P"]
                        if isinstance(self.Cell[cell]["VJ"]["P"][x], dict)
                    ],
                    index=[
                        self.Cell[cell]["VJ"]["P"][x]["sequence_id"]
                        for x in self.Cell[cell]["VJ"]["P"]
                        if isinstance(self.Cell[cell]["VJ"]["P"][x], dict)
                    ],
                )
                l_p = list(data3["sequence_id"])
                l_umi_p = [
                    int(x) for x in pd.to_numeric(data3["duplicate_count"])
                ]
                if len(l_p) > 1:
                    if "sequence_alignment" in data3:
                        l_seq_p = list(data3["sequence_alignment"])
                        if len(list(set(l_seq_p))) == 1:
                            highest_umi_l = max(l_umi_p)
                            highest_umi_l_idx = [
                                i
                                for i, j in enumerate(l_umi_p)
                                if j == highest_umi_l
                            ]
                            keep_index_l = highest_umi_l_idx[0]
                            self.drop_contig.append(
                                l_p[:keep_index_l] + l_p[keep_index_l:]
                            )
                            keep_lc_contig = l_p[keep_index_l]
                            data3.at[keep_lc_contig, "duplicate_count"] = int(
                                np.sum(
                                    l_umi_p[:keep_index_l]
                                    + l_umi_p[keep_index_l:]
                                )
                            )
                            self.umi_adjustment.update(
                                {
                                    keep_lc_contig: int(
                                        np.sum(
                                            l_umi_p[:keep_index_l]
                                            + l_umi_p[keep_index_l:]
                                        )
                                    )
                                }
                            )
                            # refresh
                            data3 = pd.DataFrame([data3.loc[keep_lc_contig]])
                            l_p = list(data3["sequence_id"])
                            l_umi_p = [
                                int(x)
                                for x in pd.to_numeric(data3["duplicate_count"])
                            ]
            if len(self.Cell[cell]["VJ"]["NP"]) > 0:
                data4 = pd.DataFrame(
                    [
                        self.Cell[cell]["VJ"]["NP"][x]
                        for x in self.Cell[cell]["VJ"]["NP"]
                        if isinstance(self.Cell[cell]["VJ"]["NP"][x], dict)
                    ],
                    index=[
                        self.Cell[cell]["VJ"]["NP"][x]["sequence_id"]
                        for x in self.Cell[cell]["VJ"]["NP"]
                        if isinstance(self.Cell[cell]["VJ"]["NP"][x], dict)
                    ],
                )
                l_np = list(data4["sequence_id"])
                l_umi_np = [
                    int(x) for x in pd.to_numeric(data4["duplicate_count"])
                ]

            if "h_p" not in locals():
                h_p = []
            if "l_p" not in locals():
                l_p = []
            if "h_np" not in locals():
                h_np = []
            if "l_np" not in locals():
                l_np = []

            if len(h_p) > 0:
                for hx in h_p:
                    v = v_dict[hx]
                    d = d_dict[hx]
                    j = j_dict[hx]
                    c = c_dict[hx]
                    if present(v):
                        if not re.search("IGH|TR[BD]|TRAV.*/DV", v):
                            self.drop_contig.append(hx)
                    if present(d):
                        if not re.search("IGH|TR[BD]", d):
                            self.drop_contig.append(hx)
                    if present(j):
                        if not re.search("IGH|TR[BD]", j):
                            self.drop_contig.append(hx)
                    if present(c):
                        if not re.search("IGH|TR[BD]", c):
                            self.drop_contig.append(hx)
                    if present(j):
                        if present(v):
                            if not_same_call(v, j, "IGH"):
                                self.drop_contig.append(hx)
                            elif not_same_call(v, j, "TRB"):
                                self.drop_contig.append(hx)
                            elif not_same_call(v, j, "TRD"):
                                if not re.search("TRAV.*/DV", v):
                                    self.drop_contig.append(hx)
                        if present(d):
                            if not_same_call(d, j, "IGH"):
                                self.drop_contig.append(hx)
                            elif not_same_call(d, j, "TRB"):
                                self.drop_contig.append(hx)
                            elif not_same_call(d, j, "TRD"):
                                self.drop_contig.append(hx)
                    else:
                        self.drop_contig.append(hx)

            if len(h_np) > 0:
                for hx in h_np:
                    v = v_dict[hx]
                    d = d_dict[hx]
                    j = j_dict[hx]
                    c = c_dict[hx]
                    if present(v):
                        if not re.search("IGH|TR[BD]|TRAV.*/DV", v):
                            self.drop_contig.append(hx)
                    if present(d):
                        if not re.search("IGH|TR[BD]", d):
                            self.drop_contig.append(hx)
                    if present(j):
                        if not re.search("IGH|TR[BD]", j):
                            self.drop_contig.append(hx)
                    if present(c):
                        if not re.search("IGH|TR[BD]", c):
                            self.drop_contig.append(hx)

                    if present(j):
                        if present(v):
                            if not_same_call(v, j, "IGH"):
                                self.drop_contig.append(hx)
                            elif not_same_call(v, j, "TRB"):
                                self.drop_contig.append(hx)
                            elif not_same_call(v, j, "TRD"):
                                if not re.search("TRAV.*/DV", v):
                                    self.drop_contig.append(hx)
                        if present(d):
                            if not_same_call(d, j, "IGH"):
                                self.drop_contig.append(hx)
                            elif not_same_call(d, j, "TRB"):
                                self.drop_contig.append(hx)
                            elif not_same_call(d, j, "TRD"):
                                self.drop_contig.append(hx)
                    else:
                        self.drop_contig.append(hx)
            if len(l_p) > 0:
                for lx in l_p:
                    v = v_dict[lx]
                    j = j_dict[lx]
                    c = c_dict[lx]
                    if present(v):
                        if re.search("IGH|TR[BD]", v):
                            self.drop_contig.append(lx)
                    if present(j):
                        if re.search("IGH|TR[BD]", j):
                            self.drop_contig.append(lx)
                    if present(c):
                        if re.search("IGH|TR[BD]", c):
                            self.drop_contig.append(lx)

                    if present(j):
                        if present(v):
                            if not_same_call(v, j, "IGK"):
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "IGL"):
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "TRA"):
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "TRG"):
                                self.drop_contig.append(lx)
                    else:
                        self.drop_contig.append(lx)

            if len(l_np) > 0:
                for lx in l_np:
                    v = v_dict[lx]
                    j = j_dict[lx]
                    c = c_dict[lx]
                    if present(v):
                        if re.search("IGH|TR[BD]", v):
                            self.drop_contig.append(lx)
                    if present(j):
                        if re.search("IGH|TR[BD]", j):
                            self.drop_contig.append(lx)
                    if present(c):
                        if re.search("IGH|TR[BD]", c):
                            self.drop_contig.append(lx)

                    if present(j):
                        if present(v):
                            if not_same_call(v, j, "IGK"):
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "IGL"):
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "TRA"):
                                self.drop_contig.append(lx)
                            elif not_same_call(v, j, "TRG"):
                                self.drop_contig.append(lx)
                    else:
                        self.drop_contig.append(lx)


def run_igblastn(
    fasta: Union[str, PathLike],
    igblast_db: Optional[str] = None,
    org: Literal["human", "mouse"] = "human",
    loci: Literal["ig", "tr"] = "ig",
    evalue: float = 1e-4,
    min_d_match: int = 9,
    verbose: bool = False,
):
    """
    Reannotate with IgBLASTn.

    Parameters
    ----------
    fasta : PathLike
        fasta file for reannotation.
    igblast_db : PathLike, Optional
        path to igblast database.
    org : str
        organism for germline sequences.
    loci : str
        `ig` or `tr` mode for running igblastn.
    evalue : float
        This is the statistical significance threshold for reporting matches
        against database sequences. Lower EXPECT thresholds are more stringent
        and report only high similarity matches. Choose higher EXPECT value
        (for example 1 or more) if you expect a low identity between your query
        sequence and the targets.
    min_d_match : int
        minimum D nucleotide match.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.

    """
    env = os.environ.copy()
    if igblast_db is None:
        try:
            igdb = env["IGDATA"]
        except KeyError:
            raise KeyError(
                (
                    "Environmental variable IGDATA must be set. Otherwise,"
                    + " please provide path to igblast database"
                )
            )
    else:
        env["IGDATA"] = igblast_db
        igdb = env["IGDATA"]

    outfolder = os.path.abspath(os.path.dirname(fasta)) + "/tmp"
    informat_dict = {"blast": "_igblast.fmt7", "airr": "_igblast.tsv"}
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    loci_type = {"ig": "Ig", "tr": "TCR"}
    outformat = {"blast": "7 std qseq sseq btop", "airr": "19"}

    for fileformat in ["blast", "airr"]:
        outfile = (
            os.path.basename(fasta).split(".fasta")[0]
            + informat_dict[fileformat]
        )
        if loci == "tr":
            cmd = [
                "igblastn",
                "-germline_db_V",
                igdb + "/database/imgt_" + org + "_" + loci + "_v",
                "-germline_db_D",
                igdb + "/database/imgt_" + org + "_" + loci + "_d",
                "-germline_db_J",
                igdb + "/database/imgt_" + org + "_" + loci + "_j",
                "-auxiliary_data",
                igdb + "optional_file/" + org + "_gl.aux",
                "-domain_system",
                "imgt",
                "-ig_seqtype",
                loci_type[loci],
                "-organism",
                org,
                "-outfmt",
                outformat[fileformat],
                "-query",
                fasta,
                "-out",
                "{}/{}".format(outfolder, outfile),
                "-evalue",
                str(evalue),
                "-min_D_match",
                str(min_d_match),
                "-D_penalty",
                str(-4),
            ]
        else:
            cmd = [
                "igblastn",
                "-germline_db_V",
                igdb + "/database/imgt_" + org + "_" + loci + "_v",
                "-germline_db_D",
                igdb + "/database/imgt_" + org + "_" + loci + "_d",
                "-germline_db_J",
                igdb + "/database/imgt_" + org + "_" + loci + "_j",
                "-auxiliary_data",
                igdb + "optional_file/" + org + "_gl.aux",
                "-domain_system",
                "imgt",
                "-ig_seqtype",
                loci_type[loci],
                "-organism",
                org,
                "-outfmt",
                outformat[fileformat],
                "-query",
                fasta,
                "-out",
                "{}/{}".format(outfolder, outfile),
                "-evalue",
                str(evalue),
                "-min_D_match",
                str(min_d_match),
            ]

        if verbose:
            print("Running command: %s\n" % (" ".join(cmd)))
        run(cmd, env=env)  # logs are printed to terminal


def assign_DJ(
    fasta: Union[str, PathLike],
    org: Literal["human", "mouse"] = "human",
    loci: Literal["ig", "tr"] = "tr",
    call: Literal["d", "j"] = "j",
    database: Optional[str] = None,
    evalue: float = 1e-4,
    max_hsps: int = 10,
    dust: Optional[Union[Literal["yes", "no"], str]] = None,
    word_size: Optional[int] = None,
    outfmt: str = (
        "6 qseqid sseqid pident length mismatch gapopen "
        + "qstart qend sstart send evalue bitscore qseq sseq"
    ),
    filename_prefix: Optional[str] = None,
    overwrite: bool = False,
    verbose: bool = False,
):
    """
    Annotate contigs with constant region call using blastn.

    Parameters
    ----------
    fasta : str, PathLike
        path to fasta file.
    org : str
        organism of reference folder. Default is 'human'.
    loci : str
        locus. 'ig' or 'tr',
    call : str
        Either 'd' of 'j' gene.
    database : str, PathLike, Optional
        path to database.
        Defaults to `$IGDATA` environmental variable if v/d/j_call.
        Defaults to `$BLASTDB` environmental variable if c_call.
    evalue : float
        This is the statistical significance threshold for reporting matches
        against database sequences. Lower EXPECT thresholds are more stringent
        and report only high similarity matches. Choose higher EXPECT value
        (for example 1 or more) if you expect a low identity between your query
        sequence and the targets.
    max_hsps: int
        Maximum number of HSPs (alignments) to keep for any single query-subject pair.
        The HSPs shown will be the best as judged by expect value. This number should
        be an integer that is one or greater. Setting it to one will show only the best
        HSP for every query-subject pair. Only affects the output file in the tmp folder.
    dust: str
        dustmasker options. Filter query sequence with DUST
        Format: 'yes', or 'no' to disable. Accepts str.
        If None, defaults to `20 64 1`.
    word_size : int
        Word size for wordfinder algorithm (length of best perfect match).
        Must be >=4. None defaults to 4.
    filename_prefix : str, Optional
        prefix of file name preceding '_contig'. None defaults to 'filtered'.
    overwrite: bool
        whether or not to overwrite the assignments.
    verbose : bool
        whether or not to print the blast command in terminal.
        Default is False.

    Returns
    -------
        D/J genes assigned.
    """
    # main function from here
    filePath = check_filepath(
        fasta, filename_prefix=filename_prefix, endswith=".fasta"
    )
    if filePath is None:
        raise FileNotFoundError(
            (
                "Path to fasta file is unknown. Please specify "
                + "path to fasta file or folder containing fasta file."
            )
        )

    # run blast
    blast_out = run_blastn(
        fasta=filePath,
        database=database,
        org=org,
        loci=loci,
        call=call,
        max_hsps=max_hsps,
        evalue=evalue,
        outfmt=outfmt,
        dust=dust,
        word_size=word_size,
        verbose=verbose,
    )

    # read the original object
    passfile = "{}/tmp/{}.tsv".format(
        os.path.dirname(filePath),
        os.path.basename(filePath).split(".fasta")[0] + "_igblast_db-pass",
    )
    failfile = "{}/tmp/{}.tsv".format(
        os.path.dirname(filePath),
        os.path.basename(filePath).split(".fasta")[0] + "_igblast_db-fail",
    )

    transfer_assignment(
        passfile=passfile,
        failfile=failfile,
        blast_result=blast_out.drop_duplicates(
            subset="sequence_id", keep="first"
        ),
        eval_threshold=evalue,
        call=call,
        overwrite=overwrite,
    )


def run_blastn(
    fasta: Union[PathLike, str],
    database: Optional[Union[PathLike, str]],
    org: Literal["human", "mouse"] = "human",
    loci: Literal["ig", "tr"] = "ig",
    call: Literal["v", "d", "j", "c"] = "c",
    max_hsps: int = 10,
    evalue: float = 1e-4,
    outfmt: str = (
        "6 qseqid sseqid pident length mismatch gapopen "
        + "qstart qend sstart send evalue bitscore qseq sseq"
    ),
    dust: Optional[Union[Literal["yes", "no"], str]] = None,
    word_size: Optional[int] = None,
    verbose: bool = False,
):
    """
    Annotate contigs using blastn.

    Parameters
    ----------
    fasta : str, PathLike
        path to fasta file.
    database : str, PathLike, Optional
        path to database.
        Defaults to `$IGDATA` environmental variable if v/d/j_call.
        Defaults to `$BLASTDB` environmental variable if c_call.
    org : str
        organism of reference folder. Default is 'human'.
    loci : str
        locus. 'ig' or 'tr',
    call : str
        Either 'v', 'd', 'j' or 'c' gene.
    max_hsps: int
        Maximum number of HSPs (alignments) to keep for any single query-subject pair.
        The HSPs shown will be the best as judged by expect value. This number should
        be an integer that is one or greater. Setting it to one will show only the best
        HSP for every query-subject pair. Only affects the output file in the tmp folder.
    evalue : float
        This is the statistical significance threshold for reporting matches
        against database sequences. Lower EXPECT thresholds are more stringent
        and report only high similarity matches. Choose higher EXPECT value
        (for example 1 or more) if you expect a low identity between your query
        sequence and the targets.
    outfmt : str
        blastn output format.
    dust: str
        dustmasker options. Filter query sequence with DUST
        Format: 'yes', or 'no' to disable. Accepts str.
        If None, defaults to `20 64 1`.
    word_size : int
        Word size for wordfinder algorithm (length of best perfect match).
        Must be >=4. None defaults to 4.
    verbose : bool
        whether or not to print the blast command in terminal.
        Default is False.

    Returns
    -------
        blastn assignment.
    """
    env = os.environ.copy()
    if call != "c":
        if database is None:
            try:
                bdb = env["IGDATA"]
            except KeyError:
                raise KeyError(
                    (
                        "Environmental variable IGDATA must be set. "
                        + "Otherwise, please provide path to igblast database."
                    )
                )
            bdb = bdb + "database/imgt_" + org + "_" + loci + "_" + call
        else:
            env["IGDATA"] = database
            bdb = database
            if not bdb.endswith("_" + loci + "_" + call):
                bdb = bdb + "database/imgt_" + org + "_" + loci + "_" + call
    else:
        if database is None:
            try:
                bdb = env["BLASTDB"]
            except KeyError:
                raise KeyError(
                    (
                        "Environmental variable BLASTDB must be set. "
                        + "Otherwise, please provide path to blast database"
                    )
                )
            bdb = bdb + org + "/" + org + "_BCR_C.fasta"
        else:
            env["BLASTDB"] = database
            bdb = database

    cmd = [
        "blastn",
        "-db",
        bdb,
        "-evalue",
        str(evalue),
        "-max_hsps",
        str(max_hsps),
        "-outfmt",
        outfmt,
        "-query",
        fasta,
    ]

    if dust is not None:
        cmd = cmd + ["-dust", str(dust)]
    if word_size is not None:
        cmd = cmd + ["-word_size", str(word_size)]

    blast_out = "{}/tmp/{}.tsv".format(
        os.path.dirname(fasta),
        os.path.basename(fasta).split(".fasta")[0] + "_" + call + "_blast",
    )

    if verbose:
        print("Running command: %s\n" % (" ".join(cmd)))
    with open(blast_out, "w") as out:
        run(cmd, stdout=out, env=env)
    try:
        dat = pd.read_csv(blast_out, sep="\t", header=None)
        dat.columns = [
            "sequence_id",
            call + "_call",
            call + "_identity",
            call + "_alignment_length",
            call + "_number_of_mismatches",
            call + "_number_of_gap_openings",
            call + "_sequence_start",
            call + "_sequence_end",
            call + "_germline_start",
            call + "_germline_end",
            call + "_support",
            call + "_score",
            call + "_sequence_alignment",
            call + "_germline_alignment",
        ]
    except pd.errors.EmptyDataError:
        dat = pd.DataFrame(
            columns=[
                "sequence_id",
                call + "_call",
                call + "_identity",
                call + "_alignment_length",
                call + "_number_of_mismatches",
                call + "_number_of_gap_openings",
                call + "_sequence_start",
                call + "_sequence_end",
                call + "_germline_start",
                call + "_germline_end",
                call + "_support",
                call + "_score",
                call + "_sequence_alignment",
                call + "_germline_alignment",
            ]
        )
    write_blastn(dat, blast_out)
    dat = load_data(dat)
    return dat


def transfer_assignment(
    passfile: Union[PathLike, str],
    failfile: Union[PathLike, str],
    blast_result: pd.DataFrame,
    eval_threshold: float,
    call: Literal["v", "d", "j", "c"] = "c",
    overwrite: bool = False,
):
    """Update gene calls with blastn results."""
    if os.path.isfile(passfile):
        db_pass = load_data(passfile)
    else:
        db_pass = None
    if os.path.isfile(failfile):
        db_fail = load_data(failfile)
        # should be pretty safe to fill this in
        db_fail["vj_in_frame"].fillna(value="F", inplace=True)
        db_fail["productive"].fillna(value="F", inplace=True)
        db_fail["c_call"].fillna(value="", inplace=True)
        db_fail["v_call"].fillna(value="", inplace=True)
        db_fail["d_call"].fillna(value="", inplace=True)
        db_fail["j_call"].fillna(value="", inplace=True)
        db_fail["locus"].fillna(value="", inplace=True)
        for i, r in db_fail.iterrows():
            if not present(r.locus):
                calls = list(
                    set(
                        [r.v_call[:3], r.d_call[:3], r.j_call[:3], r.c_call[:3]]
                    )
                )
                locus = "".join([c for c in calls if present(c)])
                if len(locus) == 3:
                    db_fail.at[i, "locus"] = locus
    else:
        db_fail = None
    if blast_result.shape[0] < 1:
        blast_result = None

    if blast_result is not None:
        if db_pass is not None:
            if call + "_support" in db_pass:
                db_pass_evalues = dict(db_pass[call + "_support"])
            if call + "_score" in db_pass:
                db_pass_scores = dict(db_pass[call + "_score"])
            db_pass[call + "_call"].fillna(value="", inplace=True)
            db_pass_call = dict(db_pass[call + "_call"])
            if call + "_support" in db_pass:
                db_pass[call + "_support_igblastn"] = pd.Series(db_pass_evalues)
            if call + "_score" in db_pass:
                db_pass[call + "_score_igblastn"] = pd.Series(db_pass_scores)
            db_pass[call + "_call_igblastn"] = pd.Series(db_pass_call)
            db_pass[call + "_call_igblastn"].fillna(value="", inplace=True)
            for col in blast_result:
                if col != "sequence_id":
                    db_pass[col + "_blastn"] = pd.Series(blast_result[col])
                    if col in [
                        call + "_call",
                        call + "_sequence_alignment",
                        call + "_germline_alignment",
                    ]:
                        db_pass[col + "_blastn"].fillna(value="", inplace=True)
            db_pass[call + "_source"] = ""
            if overwrite:
                for i in db_pass["sequence_id"]:
                    vend = db_pass.loc[i, "v_sequence_end"]
                    if not present(vend):
                        vend_ = 0
                    else:
                        vend_ = vend
                    jstart = db_pass.loc[i, "j_sequence_start"]
                    if not present(jstart):
                        jstart_ = 1000
                    else:
                        jstart_ = jstart
                    callstart = db_pass.loc[i, call + "_sequence_start_blastn"]
                    callend = db_pass.loc[i, call + "_sequence_end_blastn"]
                    if (callstart >= vend_) and (callend <= jstart_):
                        if call + "_support_igblastn" in db_pass:
                            eval1 = db_pass.loc[i, call + "_support_igblastn"]
                        else:
                            eval1 = 1
                        eval2 = db_pass.loc[i, call + "_support_blastn"]
                        if (
                            db_pass.loc[i, call + "_call_igblastn"]
                            != db_pass.loc[i, call + "_call_blastn"]
                        ):
                            if call + "_call_10x" in db_pass:
                                if (
                                    re.sub(
                                        "[*][0-9][0-9]",
                                        "",
                                        db_pass.loc[i, call + "_call_blastn"],
                                    )
                                    != db_pass.loc[i, call + "_call_10x"]
                                ):
                                    if present(eval1):
                                        if eval1 > eval2:
                                            db_pass.at[
                                                i, call + "_call"
                                            ] = db_pass.at[
                                                i, call + "_call_blastn"
                                            ]
                                            db_pass.at[
                                                i, call + "_sequence_start"
                                            ] = db_pass.at[
                                                i,
                                                call + "_sequence_start_blastn",
                                            ]
                                            db_pass.at[
                                                i, call + "_sequence_end"
                                            ] = db_pass.at[
                                                i, call + "_sequence_end_blastn"
                                            ]
                                            db_pass.at[
                                                i, call + "_germline_start"
                                            ] = db_pass.at[
                                                i,
                                                call + "_germline_start_blastn",
                                            ]
                                            db_pass.at[
                                                i, call + "_germline_end"
                                            ] = db_pass.at[
                                                i, call + "_germline_end_blastn"
                                            ]
                                            db_pass.at[
                                                i, call + "_source"
                                            ] = "blastn"
                                    else:
                                        if present(eval2):
                                            db_pass.at[
                                                i, call + "_call"
                                            ] = db_pass.at[
                                                i, call + "_call_blastn"
                                            ]
                                            db_pass.at[
                                                i, call + "_sequence_start"
                                            ] = db_pass.at[
                                                i,
                                                call + "_sequence_start_blastn",
                                            ]
                                            db_pass.at[
                                                i, call + "_sequence_end"
                                            ] = db_pass.at[
                                                i, call + "_sequence_end_blastn"
                                            ]
                                            db_pass.at[
                                                i, call + "_germline_start"
                                            ] = db_pass.at[
                                                i,
                                                call + "_germline_start_blastn",
                                            ]
                                            db_pass.at[
                                                i, call + "_germline_end"
                                            ] = db_pass.at[
                                                i, call + "_germline_end_blastn"
                                            ]
                                            db_pass.at[
                                                i, call + "_source"
                                            ] = "blastn"
                                else:
                                    db_pass.at[i, call + "_source"] = "10x"
                                    db_pass.at[i, call + "_call"] = db_pass.at[
                                        i, call + "_call_blastn"
                                    ]
                                    if present(db_pass.loc[i, "junction_10x"]):
                                        if present(db_pass.loc[i, "junction"]):
                                            if (
                                                db_pass.loc[i, "junction"]
                                                != db_pass.loc[
                                                    i, "junction_10x"
                                                ]
                                            ):
                                                db_pass.at[
                                                    i, "junction"
                                                ] = db_pass.at[
                                                    i, "junction_10x"
                                                ]
                                                db_pass.at[
                                                    i, "junction_aa"
                                                ] = db_pass.at[
                                                    i, "junction_10x_aa"
                                                ]
                        else:
                            if present(eval1):
                                if eval1 > eval2:
                                    db_pass.at[i, call + "_call"] = db_pass.at[
                                        i, call + "_call_blastn"
                                    ]
                                    db_pass.at[
                                        i, call + "_sequence_start"
                                    ] = db_pass.at[
                                        i, call + "_sequence_start_blastn"
                                    ]
                                    db_pass.at[
                                        i, call + "_sequence_end"
                                    ] = db_pass.at[
                                        i, call + "_sequence_end_blastn"
                                    ]
                                    db_pass.at[
                                        i, call + "_germline_start"
                                    ] = db_pass.at[
                                        i, call + "_germline_start_blastn"
                                    ]
                                    db_pass.at[
                                        i, call + "_germline_end"
                                    ] = db_pass.at[
                                        i, call + "_germline_end_blastn"
                                    ]
                                    db_pass.at[i, call + "_source"] = "blastn"
                            else:
                                if present(eval2):
                                    db_pass.at[i, call + "_call"] = db_pass.at[
                                        i, call + "_call_blastn"
                                    ]
                                    db_pass.at[
                                        i, call + "_sequence_start"
                                    ] = db_pass.at[
                                        i, call + "_sequence_start_blastn"
                                    ]
                                    db_pass.at[
                                        i, call + "_sequence_end"
                                    ] = db_pass.at[
                                        i, call + "_sequence_end_blastn"
                                    ]
                                    db_pass.at[
                                        i, call + "_germline_start"
                                    ] = db_pass.at[
                                        i, call + "_germline_start_blastn"
                                    ]
                                    db_pass.at[
                                        i, call + "_germline_end"
                                    ] = db_pass.at[
                                        i, call + "_germline_end_blastn"
                                    ]
                                    db_pass.at[i, call + "_source"] = "blastn"

                vend = db_pass["v_sequence_end"]
                dstart = db_pass["d_sequence_start"]
                dend = db_pass["d_sequence_end"]
                jstart = db_pass["j_sequence_start"]

                np1 = [
                    str(int(n)) if n >= 0 else ""
                    for n in [
                        (d - v) - 1
                        if pd.notnull(v) and pd.notnull(d)
                        else np.nan
                        for v, d in zip(vend, dstart)
                    ]
                ]
                np2 = [
                    str(int(n)) if n >= 0 else ""
                    for n in [
                        (j - d) - 1
                        if pd.notnull(j) and pd.notnull(d)
                        else np.nan
                        for d, j in zip(dend, jstart)
                    ]
                ]

                db_pass["np1_length"] = np1
                db_pass["np2_length"] = np2

                for i in db_pass["sequence_id"]:
                    if not present(db_pass.loc[i, "np1_length"]):
                        vend = db_pass.loc[i, "v_sequence_end"]
                        if present(vend):
                            jstart = db_pass.loc[i, "j_sequence_start"]
                            if present(jstart):
                                np1l = (jstart - vend) - 1
                                if np1l >= 0:
                                    db_pass.loc[i, "np1_length"] = np1l
            # fill in blanks
            db_pass = sanitize_data(db_pass)
            db_pass.to_csv(passfile, sep="\t", index=False)

        if db_fail is not None:
            if call + "_support" in db_fail:
                db_fail_evalues = dict(db_fail[call + "_support"])
            if call + "_score" in db_fail:
                db_fail_scores = dict(db_fail[call + "_score"])
            db_fail[call + "_call"].fillna(value="", inplace=True)
            db_fail_call = dict(db_fail[call + "_call"])
            if call + "_support" in db_fail:
                db_fail[call + "_support_igblastn"] = pd.Series(db_fail_evalues)
            if call + "_score" in db_fail:
                db_fail[call + "_score_igblastn"] = pd.Series(db_fail_scores)
            db_fail[call + "_call_igblastn"] = pd.Series(db_fail_call)
            db_fail[call + "_call_igblastn"].fillna(value="", inplace=True)
            for col in blast_result:
                if col != "sequence_id":
                    db_fail[col + "_blastn"] = pd.Series(blast_result[col])
                    if col in [
                        call + "_call",
                        call + "_sequence_alignment",
                        call + "_germline_alignment",
                    ]:
                        db_fail[col + "_blastn"].fillna(value="", inplace=True)
            db_fail[call + "_source"] = ""
            if overwrite:
                for i in db_fail["sequence_id"]:
                    vend = db_fail.loc[i, "v_sequence_end"]
                    if not present(vend):
                        vend_ = 0
                    else:
                        vend_ = vend
                    jstart = db_fail.loc[i, "j_sequence_start"]
                    if not present(jstart):
                        jstart_ = 1000
                    else:
                        jstart_ = jstart
                    callstart = db_fail.loc[i, call + "_sequence_start_blastn"]
                    callend = db_fail.loc[i, call + "_sequence_end_blastn"]
                    if (callstart >= vend_) and (callend <= jstart_):
                        if call + "_support_igblastn" in db_fail:
                            eval1 = db_fail.loc[i, call + "_support_igblastn"]
                        else:
                            eval1 = 1
                        eval2 = db_fail.loc[i, call + "_support_blastn"]
                        if (
                            db_fail.loc[i, call + "_call_igblastn"]
                            != db_fail.loc[i, call + "_call_blastn"]
                        ):
                            if call + "_call_10x" in db_fail:
                                if (
                                    re.sub(
                                        "[*][0-9][0-9]",
                                        "",
                                        db_fail.loc[i, call + "_call_blastn"],
                                    )
                                    != db_fail.loc[i, call + "_call_10x"]
                                ):
                                    if present(eval1):
                                        if eval1 > eval2:
                                            db_fail.at[
                                                i, call + "_call"
                                            ] = db_fail.at[
                                                i, call + "_call_blastn"
                                            ]
                                            db_fail.at[
                                                i, call + "_sequence_start"
                                            ] = db_fail.at[
                                                i,
                                                call + "_sequence_start_blastn",
                                            ]
                                            db_fail.at[
                                                i, call + "_sequence_end"
                                            ] = db_fail.at[
                                                i, call + "_sequence_end_blastn"
                                            ]
                                            db_fail.at[
                                                i, call + "_germline_start"
                                            ] = db_fail.at[
                                                i,
                                                call + "_germline_start_blastn",
                                            ]
                                            db_fail.at[
                                                i, call + "_germline_end"
                                            ] = db_fail.at[
                                                i, call + "_germline_end_blastn"
                                            ]
                                            db_fail.at[
                                                i, call + "_source"
                                            ] = "blastn"
                                    else:
                                        if present(eval2):
                                            db_fail.at[
                                                i, call + "_call"
                                            ] = db_fail.at[
                                                i, call + "_call_blastn"
                                            ]
                                            db_fail.at[
                                                i, call + "_sequence_start"
                                            ] = db_fail.at[
                                                i,
                                                call + "_sequence_start_blastn",
                                            ]
                                            db_fail.at[
                                                i, call + "_sequence_end"
                                            ] = db_fail.at[
                                                i, call + "_sequence_end_blastn"
                                            ]
                                            db_fail.at[
                                                i, call + "_germline_start"
                                            ] = db_fail.at[
                                                i,
                                                call + "_germline_start_blastn",
                                            ]
                                            db_fail.at[
                                                i, call + "_germline_end"
                                            ] = db_fail.at[
                                                i, call + "_germline_end_blastn"
                                            ]
                                            db_fail.at[
                                                i, call + "_source"
                                            ] = "blastn"
                                else:
                                    db_fail.at[i, call + "_source"] = "10x"
                                    db_fail.at[i, call + "_call"] = db_fail.at[
                                        i, call + "_call_blastn"
                                    ]
                                    if present(db_fail.loc[i, "junction_10x"]):
                                        if present(db_fail.loc[i, "junction"]):
                                            if (
                                                db_fail.loc[i, "junction"]
                                                != db_fail.loc[
                                                    i, "junction_10x"
                                                ]
                                            ):
                                                db_fail.at[
                                                    i, "junction"
                                                ] = db_fail.at[
                                                    i, "junction_10x"
                                                ]
                                                db_fail.at[
                                                    i, "junction_aa"
                                                ] = db_fail.at[
                                                    i, "junction_10x_aa"
                                                ]
                        else:
                            if present(eval1):
                                if eval1 > eval2:
                                    db_fail.at[i, call + "_call"] = db_fail.at[
                                        i, call + "_call_blastn"
                                    ]
                                    db_fail.at[
                                        i, call + "_sequence_start"
                                    ] = db_fail.at[
                                        i, call + "_sequence_start_blastn"
                                    ]
                                    db_fail.at[
                                        i, call + "_sequence_end"
                                    ] = db_fail.at[
                                        i, call + "_sequence_end_blastn"
                                    ]
                                    db_fail.at[
                                        i, call + "_germline_start"
                                    ] = db_fail.at[
                                        i, call + "_germline_start_blastn"
                                    ]
                                    db_fail.at[
                                        i, call + "_germline_end"
                                    ] = db_fail.at[
                                        i, call + "_germline_end_blastn"
                                    ]
                                    db_fail.at[i, call + "_source"] = "blastn"
                            else:
                                if present(eval2):
                                    db_fail.at[i, call + "_call"] = db_fail.at[
                                        i, call + "_call_blastn"
                                    ]
                                    db_fail.at[
                                        i, call + "_sequence_start"
                                    ] = db_fail.at[
                                        i, call + "_sequence_start_blastn"
                                    ]
                                    db_fail.at[
                                        i, call + "_sequence_end"
                                    ] = db_fail.at[
                                        i, call + "_sequence_end_blastn"
                                    ]
                                    db_fail.at[
                                        i, call + "_germline_start"
                                    ] = db_fail.at[
                                        i, call + "_germline_start_blastn"
                                    ]
                                    db_fail.at[
                                        i, call + "_germline_end"
                                    ] = db_fail.at[
                                        i, call + "_germline_end_blastn"
                                    ]
                                    db_fail.at[i, call + "_source"] = "blastn"

                vend = db_fail["v_sequence_end"]
                dstart = db_fail["d_sequence_start"]
                dend = db_fail["d_sequence_end"]
                jstart = db_fail["j_sequence_start"]

                np1 = [
                    str(int(n)) if n >= 0 else ""
                    for n in [
                        (d - v) - 1
                        if pd.notnull(v) and pd.notnull(d)
                        else np.nan
                        for v, d in zip(vend, dstart)
                    ]
                ]
                np2 = [
                    str(int(n)) if n >= 0 else ""
                    for n in [
                        (j - d) - 1
                        if pd.notnull(j) and pd.notnull(d)
                        else np.nan
                        for d, j in zip(dend, jstart)
                    ]
                ]
                db_fail["np1_length"] = np1
                db_fail["np2_length"] = np2

                # rescue the d blanks
                for i in db_fail["sequence_id"]:
                    if not present(db_fail.loc[i, "np1_length"]):
                        vend = db_fail.loc[i, "v_sequence_end"]
                        if present(vend):
                            jstart = db_fail.loc[i, "j_sequence_start"]
                            if present(jstart):
                                np1l = (jstart - vend) - 1
                                if np1l >= 0:
                                    db_fail.loc[i, "np1_length"] = np1l

            # fill in blanks
            db_fail = sanitize_data(db_fail)
            db_fail.to_csv(failfile, sep="\t", index=False)
