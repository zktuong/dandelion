#!/usr/bin/env python
import pandas as pd
import dandelion as ddl
from pathlib import Path

from fixtures_mouse import (fasta_10x_mouse, annotation_10x_mouse,
                            database_paths_mouse, dummy_adata_mouse,
                            balbc_ighg_primers)
from fixtures import (processed_files, create_testfolder)


def test_write_fasta(create_testfolder, fasta_10x_mouse):
    out_fasta = str(create_testfolder) + "/filtered_contig.fasta"
    fh = open(out_fasta, "w")
    fh.close()
    out = ''
    for l in fasta_10x_mouse:
        out = '>' + l + '\n' + fasta_10x_mouse[l] + '\n'
        ddl.utl.Write_output(out, out_fasta)
    assert len(list(create_testfolder.iterdir())) == 1


def test_write_annotation(create_testfolder, annotation_10x_mouse):
    out_file = str(create_testfolder) + "/filtered_contig_annotations.csv"
    annotation_10x_mouse.to_csv(out_file, index=False)
    assert len(list(create_testfolder.iterdir())) == 2


def test_formatfasta(create_testfolder):
    ddl.pp.format_fastas(str(create_testfolder))
    assert len(list((create_testfolder / 'dandelion').iterdir())) == 2


def test_reannotategenes_original(create_testfolder, database_paths_mouse):
    ddl.pp.format_fastas(str(create_testfolder))
    ddl.pp.reannotate_genes(str(create_testfolder),
                            igblast_db=database_paths_mouse['igblast_db'],
                            germline=database_paths_mouse['germline'],
                            flavour = 'original',
                            org='mouse')
    assert len(list((create_testfolder / 'dandelion/tmp').iterdir())) == 4


def test_reannotategenes_other(create_testfolder, database_paths_mouse):
    ddl.pp.format_fastas(str(create_testfolder))
    ddl.pp.reannotate_genes(str(create_testfolder),
                            igblast_db=database_paths_mouse['igblast_db'],
                            germline=database_paths_mouse['germline'],
                            extended = False,
                            org='mouse')
    assert len(list((create_testfolder / 'dandelion/tmp').iterdir())) == 6
