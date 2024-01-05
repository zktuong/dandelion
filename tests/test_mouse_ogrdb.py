#!/usr/bin/env python
import pandas as pd
import dandelion as ddl
import pytest


@pytest.mark.usefixtures("create_testfolder", "fasta_10x_mouse")
def test_write_fasta(create_testfolder, fasta_10x_mouse):
    """test_write_fasta"""
    out_fasta = create_testfolder / "filtered_contig.fasta"
    ddl.utl.write_fasta(fasta_dict=fasta_10x_mouse, out_fasta=out_fasta)
    assert len(list(create_testfolder.iterdir())) == 1


@pytest.mark.usefixtures("create_testfolder", "annotation_10x_mouse")
def test_write_annotation(create_testfolder, annotation_10x_mouse):
    """test_write_annotation"""
    out_file = create_testfolder / "filtered_contig_annotations.csv"
    annotation_10x_mouse.to_csv(out_file, index=False)
    assert len(list(create_testfolder.iterdir())) == 2


@pytest.mark.usefixtures("create_testfolder")
def test_formatfasta(create_testfolder):
    """test_formatfasta"""
    ddl.pp.format_fastas(create_testfolder)
    assert len(list((create_testfolder / "dandelion").iterdir())) == 2


@pytest.mark.parametrize("flavour", ["strict", "original"])
@pytest.mark.usefixtures("create_testfolder", "database_paths_mouse")
def test_reannotategenes(create_testfolder, database_paths_mouse, flavour):
    """test_reannotategenes"""
    # different igblast versions/references may give different results and lead to regression.
    # disabling those checks would work for now.
    ddl.pp.reannotate_genes(
        create_testfolder,
        igblast_db=database_paths_mouse["igblast_db"],
        germline=database_paths_mouse["ogrdb"],
        org="mouse",
        db="ogrdb",
        strain=None,
        flavour=flavour,
    )


@pytest.mark.parametrize("flavour", ["strict", "original"])
@pytest.mark.usefixtures("create_testfolder", "database_paths_mouse")
def test_reannotategenes_nod(create_testfolder, database_paths_mouse, flavour):
    """test_reannotategenes"""
    # different igblast versions/references may give different results and lead to regression.
    # disabling those checks would work for now.
    ddl.pp.reannotate_genes(
        create_testfolder,
        igblast_db=database_paths_mouse["igblast_db"],
        germline=database_paths_mouse["ogrdb"],
        org="mouse",
        db="ogrdb",
        strain="C57BL_6J",
        flavour=flavour,
    )
