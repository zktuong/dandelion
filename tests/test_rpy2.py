#!/usr/bin/env python
import pytest

# import sys
import pandas as pd
import dandelion as ddl


# @pytest.mark.skipif(sys.platform == "darwin", reason="macos CI stalls.")
def test_importrpy2():
    """test_importrpy2"""

    from rpy2.robjects.packages import importr

    sh = importr("shazam")

    assert sh.__module__ == "rpy2.robjects.packages"


@pytest.mark.usefixtures("create_testfolder", "airr_reannotated")
# @pytest.mark.skipif(sys.platform == "darwin", reason="macos CI stalls.")
def test_mutation(create_testfolder, airr_reannotated):
    """test_mutation"""
    f = create_testfolder / "test.tsv"
    airr_reannotated.to_csv(f, sep="\t", index=False)
    ddl.pp.quantify_mutations(f)
    out = pd.read_csv(f, sep="\t")
    vdj = ddl.Dandelion(out)
    assert not vdj.data.mu_count.empty
    ddl.pp.quantify_mutations(f, frequency=True)
    assert not vdj.data.mu_freq.empty


@pytest.mark.usefixtures("create_testfolder", "database_paths")
# @pytest.mark.skipif(sys.platform == "darwin", reason="macos CI stalls.")
def test_create_germlines(create_testfolder, database_paths):
    """test create germlines"""
    f = create_testfolder / "test.tsv"
    out = pd.read_csv(f, sep="\t")
    vdj = ddl.Dandelion(out)
    ddl.pp.create_germlines(vdj, germline=database_paths["germline"])
    assert not vdj.data.germline_alignment_d_mask.empty


@pytest.mark.usefixtures("create_testfolder")
# @pytest.mark.skipif(sys.platform == "darwin", reason="macos CI stalls.")
def test_manual_threshold_and_define_clones(create_testfolder):
    """test threshold"""
    f = create_testfolder / "test.tsv"
    out = pd.read_csv(f, sep="\t")
    vdj = ddl.Dandelion(out)
    vdj.threshold = 0.1
    ddl.tl.define_clones(vdj)
    assert not vdj.data.clone_id.empty
    ddl.tl.define_clones(vdj, key_added="changeo_clone")
    assert not vdj.data.changeo_clone.empty


@pytest.mark.usefixtures("create_testfolder")
# @pytest.mark.skipif(sys.platform == "darwin", reason="macos CI stalls.")
def test_define_clones_outdir(create_testfolder):
    """test threshold"""
    f = create_testfolder / "test.tsv"
    out = pd.read_csv(f, sep="\t")
    vdj = ddl.Dandelion(out)
    vdj.threshold = 0.1
    out_path = (
        create_testfolder / "test" / "test"
    )  # this path shouldn't exist initially
    ddl.tl.define_clones(vdj, out_dir=out_path)
    assert len(list(out_path.iterdir())) == 3
    assert len(list((out_path / "tmp").iterdir())) == 2


@pytest.mark.usefixtures("create_testfolder")
# @pytest.mark.skipif(sys.platform == "darwin", reason="macos CI stalls.")
def test_scoper_i(create_testfolder):
    """test identical clones from scoper"""
    f = create_testfolder / "test.tsv"
    vdj = ddl.Dandelion(f)
    assert "clone_id" not in vdj.data
    from dandelion.external.immcantation.scoper import identical_clones

    identical_clones(vdj)
    assert not vdj.data.clone_id.empty


@pytest.mark.usefixtures("create_testfolder")
# @pytest.mark.skipif(sys.platform == "darwin", reason="macos CI stalls.")
def test_scoper_h(create_testfolder):
    """test hierarchical clones from scoper"""
    f = create_testfolder / "test.tsv"
    vdj = ddl.Dandelion(f)
    assert "clone_id" not in vdj.data
    from dandelion.external.immcantation.scoper import hierarchical_clones

    hierarchical_clones(vdj, threshold=0.15)
    assert not vdj.data.clone_id.empty


@pytest.mark.usefixtures("create_testfolder")
# @pytest.mark.skipif(sys.platform == "darwin", reason="macos CI stalls.")
def test_scoper_h(create_testfolder):
    """test spectral clones from scoper"""
    f = create_testfolder / "test.tsv"
    vdj = ddl.Dandelion(f)
    assert "clone_id" not in vdj.data
    from dandelion.external.immcantation.scoper import spectral_clones

    spectral_clones(vdj, method="novj")
    assert not vdj.data.clone_id.empty

    vdj = ddl.Dandelion(f)
    assert "clone_id" not in vdj.data
    spectral_clones(vdj, method="novj", threshold=0.15)
    assert not vdj.data.clone_id.empty

    vdj = ddl.Dandelion(f)
    assert "clone_id" not in vdj.data
    spectral_clones(vdj, method="vj", threshold=0.15)
    assert not vdj.data.clone_id.empty
