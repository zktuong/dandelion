#!/usr/bin/env python
"""test mouse"""
import pandas as pd
import dandelion as ddl
import pytest

from pathlib import Path


@pytest.mark.usefixtures("create_testfolder", "fasta_10x_mouse")
def test_write_fasta(create_testfolder, fasta_10x_mouse):
    """test_write_fasta"""
    out_fasta = str(create_testfolder) + "/filtered_contig.fasta"
    fh = open(out_fasta, "w")
    fh.close()
    out = ""
    for line in fasta_10x_mouse:
        out = ">" + line + "\n" + fasta_10x_mouse[line] + "\n"
        ddl.utl.Write_output(out, out_fasta)
    assert len(list(create_testfolder.iterdir())) == 1


@pytest.mark.usefixtures("create_testfolder", "annotation_10x_mouse")
def test_write_annotation(create_testfolder, annotation_10x_mouse):
    """test_write_annotation"""
    out_file = str(create_testfolder) + "/filtered_contig_annotations.csv"
    annotation_10x_mouse.to_csv(out_file, index=False)
    assert len(list(create_testfolder.iterdir())) == 2


@pytest.mark.usefixtures("create_testfolder")
def test_formatfasta(create_testfolder):
    """test_formatfasta"""
    ddl.pp.format_fastas(str(create_testfolder))
    assert len(list((create_testfolder / "dandelion").iterdir())) == 2


@pytest.mark.usefixtures("create_testfolder", "database_paths_mouse")
def test_reannotategenes(create_testfolder, database_paths_mouse):
    """test_reannotategenes"""
    ddl.pp.reannotate_genes(
        str(create_testfolder),
        igblast_db=database_paths_mouse["igblast_db"],
        germline=database_paths_mouse["germline"],
        org="mouse",
        reassign_dj=False,
    )
    assert len(list((create_testfolder / "dandelion" / "tmp").iterdir())) == 4


@pytest.mark.usefixtures("create_testfolder", "database_paths_mouse")
def test_reassignalleles(create_testfolder, database_paths_mouse):
    """test_reassignalleles"""
    ddl.pp.reassign_alleles(
        str(create_testfolder),
        combined_folder="test_mouse",
        germline=database_paths_mouse["germline"],
        org="mouse",
        novel=True,
        plot=False,
    )
    assert len(list((create_testfolder / "dandelion" / "tmp").iterdir())) == 7


@pytest.mark.usefixtures("database_paths_mouse")
def test_updateblastdb(database_paths_mouse):
    """test_updateblastdb"""
    ddl.utl.makeblastdb(database_paths_mouse["blastdb_fasta"])
    # assert len(list(Path(database_paths_mouse["blastdb"]).iterdir())) == 10


@pytest.mark.usefixtures(
    "create_testfolder", "database_paths_mouse", "balbc_ighg_primers"
)
def test_assignsisotypes(
    create_testfolder, database_paths_mouse, balbc_ighg_primers
):
    """test_assignsisotypes"""
    ddl.pp.assign_isotypes(
        str(create_testfolder),
        blastdb=database_paths_mouse["blastdb_fasta"],
        correction_dict=balbc_ighg_primers,
        plot=False,
    )
    assert len(list((create_testfolder / "dandelion").iterdir())) == 2


@pytest.mark.usefixtures(
    "create_testfolder", "processed_files", "database_paths_mouse"
)
def test_create_germlines(
    create_testfolder, processed_files, database_paths_mouse
):
    """test_create_germlines"""
    f = create_testfolder / "dandelion" / processed_files["filtered"]
    ddl.pp.create_germlines(f, germline=database_paths_mouse["germline"])
    f2 = create_testfolder / str(
        "dandelion/" + processed_files["germline-dmask"]
    )
    dat = pd.read_csv(f2, sep="\t")
    assert not dat["germline_alignment_d_mask"].empty


@pytest.mark.usefixtures(
    "create_testfolder", "processed_files", "dummy_adata_mouse"
)
def test_filtercontigs(create_testfolder, processed_files, dummy_adata_mouse):
    """test_filtercontigs"""
    f = create_testfolder / "dandelion" / processed_files["filtered"]
    dat = pd.read_csv(f, sep="\t")
    vdj, adata = ddl.pp.filter_contigs(dat, dummy_adata_mouse)
    f1 = create_testfolder / "test.h5"
    f2 = create_testfolder / "test.h5ad"
    vdj.write_h5(f1)
    adata.write_h5ad(f2)
    assert dat.shape[0] == 1285
    assert vdj.data.shape[0] == 956
    assert vdj.metadata.shape[0] == 448
    assert adata.n_obs == 547


@pytest.mark.usefixtures("create_testfolder")
def test_generate_network(create_testfolder):
    """test_generate_network"""
    f = create_testfolder / "test.h5"
    vdj = ddl.read_h5(f)
    with pytest.raises(ValueError):
        ddl.tl.generate_network(vdj, compute_layout=False)
    ddl.tl.find_clones(vdj)
    ddl.tl.generate_network(vdj, layout_method="mod_fr")
    assert vdj.n_obs == 448
    assert vdj.layout is not None
    assert vdj.graph is not None
    assert len(vdj.graph[1]) == 8
    ddl.tl.generate_network(vdj, compute_layout=False, min_size=3)
    assert len(vdj.graph[1]) == 0


@pytest.mark.usefixtures(
    "create_testfolder", "processed_files", "dummy_adata_mouse"
)
def test_filtercontigs_drop_contigs(
    create_testfolder, processed_files, dummy_adata_mouse
):
    """test_filtercontigs_drop_contigs"""
    f = create_testfolder / "dandelion" / processed_files["filtered"]
    dat = pd.read_csv(f, sep="\t")
    vdj, adata = ddl.pp.filter_contigs(
        dat, dummy_adata_mouse, filter_poorqualitycontig=True
    )
    assert dat.shape[0] == 1285
    assert vdj.data.shape[0] == 955
    assert vdj.metadata.shape[0] == 447
    assert adata.n_obs == 547
