#!/usr/bin/env python
import pandas as pd
import dandelion as ddl
import pytest


@pytest.mark.usefixtures("create_testfolder", "airr_reannotated", "dummy_adata")
def test_setup(create_testfolder, airr_reannotated, dummy_adata):
    """test setup"""
    vdj, adata = ddl.pp.check_contigs(airr_reannotated, dummy_adata)
    assert airr_reannotated.shape[0] == 8
    assert vdj.data.shape[0] == 8
    assert vdj.metadata.shape[0] == 5
    assert adata.n_obs == 5
    f = create_testfolder / "test.h5ddl"
    vdj.write_h5ddl(f)
    assert len(list(create_testfolder.iterdir())) == 1
    vdj2 = ddl.read_h5ddl(f)
    assert vdj2.metadata is not None


@pytest.mark.usefixtures("create_testfolder")
def test_find_clones(create_testfolder):
    """test find clones"""
    f = create_testfolder / "test.h5ddl"
    vdj = ddl.read_h5ddl(f)
    ddl.tl.find_clones(vdj)
    assert not vdj.data.clone_id.empty
    assert not vdj.metadata.clone_id.empty
    vdj.write_h5ddl(f)
    with pytest.raises(ValueError):
        ddl.tl.find_clones(vdj, key="random_column")


@pytest.mark.usefixtures("airr_reannotated")
def test_find_clonesfromfile(airr_reannotated):
    """test find clones"""
    vdj = ddl.tl.find_clones(airr_reannotated)
    assert not vdj.data.clone_id.empty
    assert not vdj.metadata.clone_id.empty
    vdj2 = ddl.tl.find_clones(airr_reannotated, by_alleles=True)
    assert not vdj2.data.clone_id.empty
    assert not vdj2.metadata.clone_id.empty
