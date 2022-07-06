#!/usr/bin/env python
"""test tools 2"""
import pandas as pd
import dandelion as ddl
import pytest


@pytest.mark.usefixtures("create_testfolder", "airr_reannotated", "dummy_adata")
def test_setup(create_testfolder, airr_reannotated, dummy_adata):
    """test setup"""
    vdj, adata = ddl.pp.filter_contigs(airr_reannotated, dummy_adata)
    assert airr_reannotated.shape[0] == 8
    assert vdj.data.shape[0] == 7
    assert vdj.metadata.shape[0] == 4
    assert adata.n_obs == 5
    f = create_testfolder / "test.h5"
    vdj.write_h5(f)
    assert len(list(create_testfolder.iterdir())) == 1
    vdj2 = ddl.read_h5(f)
    assert vdj2.metadata is not None


@pytest.mark.usefixtures("create_testfolder")
def test_find_clones(create_testfolder):
    """test find clones"""
    f = create_testfolder / "test.h5"
    vdj = ddl.read_h5(f)
    ddl.tl.find_clones(vdj, collapse_label=True)
    assert not vdj.data.clone_id.empty
    assert not vdj.metadata.clone_id.empty
    assert len(set(x for x in vdj.metadata["clone_id"] if pd.notnull(x))) == 4
    vdj.write_h5(f)
    with pytest.raises(ValueError):
        ddl.tl.find_clones(vdj, key="random_column")


@pytest.mark.usefixtures("airr_reannotated")
def test_find_clonesfromfile(airr_reannotated):
    """test find clones"""
    vdj = ddl.tl.find_clones(airr_reannotated, collapse_label=True)
    assert not vdj.data.clone_id.empty
    assert not vdj.metadata.clone_id.empty
    vdj2 = ddl.tl.find_clones(airr_reannotated, by_alleles=True)
    assert not vdj2.data.clone_id.empty
    assert not vdj2.metadata.clone_id.empty
