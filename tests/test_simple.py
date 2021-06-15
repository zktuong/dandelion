#!/usr/bin/env python
import pytest
import pandas as pd
import dandelion as ddl
from pathlib import Path

from fixtures import (fasta_10x, annotation_10x, airr_10x, create_testfolder,
                      database_paths)


def test_write_fasta(create_testfolder, fasta_10x):
    out_fasta = str(create_testfolder) + "/filtered_contig.fasta"
    fh = open(out_fasta, "w")
    fh.close()
    out = ''
    for l in fasta_10x:
        out = '>' + l + '\n' + fasta_10x[l] + '\n'
        ddl.utl.Write_output(out, out_fasta)
    assert len(list(create_testfolder.iterdir())) == 1


def test_write_annotation(create_testfolder, annotation_10x):
    out_file = str(create_testfolder) + "/filtered_contig_annotations.csv"
    annotation_10x.to_csv(out_file, index=False)
    assert len(list(create_testfolder.iterdir())) == 2


def test_write_airr(create_testfolder, airr_10x):
    out_file = str(create_testfolder) + "/airr_rearrangements.tsv"
    airr_10x.to_csv(out_file, sep='\t', index=False)
    assert len(list(create_testfolder.iterdir())) == 3


def test_read10xairr(create_testfolder):
    airr_file = str(create_testfolder) + "/airr_rearrangements.tsv"
    vdj = ddl.read_10x_airr(airr_file)
    assert vdj.data.shape[0] == 9
    assert vdj.data.shape[1] == 33
    assert vdj.metadata.shape[0] == 4
    assert vdj.metadata.shape[1] == 23


def test_formatfasta(create_testfolder):
    ddl.pp.format_fastas(str(create_testfolder))
    assert len(list((create_testfolder / 'dandelion').iterdir())) == 2


def test_reannotategenes(create_testfolder, database_paths):
    ddl.pp.reannotate_genes(str(create_testfolder),
                            igblast_db=database_paths['igblast_db'],
                            germline=database_paths['germline'])
    assert len(list((create_testfolder / 'dandelion/tmp').iterdir())) == 3
