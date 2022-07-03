#!/usr/bin/env python
"""test load_data"""

import dandelion as ddl
import pytest


@pytest.mark.usefixtures("airr_reannotated")
def test_load_data(airr_reannotated):
    """test load_data"""
    vdj = ddl.Dandelion(airr_reannotated)
    assert all(
        [x != y for x, y in zip(vdj.data["cell_id"], vdj.data["sequence_id"])]
    )
    cell_ids = list(vdj.data["cell_id"])
    tmp = vdj.data.drop("cell_id", axis=1)
    vdj = ddl.Dandelion(tmp)
    assert all([x == y for x, y in zip(vdj.data["cell_id"], cell_ids)])


@pytest.mark.usefixtures("airr_generic")
def test_slice_data(airr_generic):
    """test load_data"""
    vdj = ddl.Dandelion(airr_generic)
    assert vdj.data.shape[0] == 105
    assert vdj.metadata.shape[0] == 40
    vdj2 = vdj[vdj.data["productive"] == "T"]
    assert vdj2.data.shape[0] == 94
    assert vdj2.metadata.shape[0] == 38
    vdj2 = vdj[vdj.metadata["productive_VDJ"] == "T"]
    assert vdj2.data.shape[0] == 36
    assert vdj2.metadata.shape[0] == 17


@pytest.mark.usefixtures("airr_generic")
def test_names(airr_generic):
    """test load_data"""
    vdj = ddl.Dandelion(airr_generic)
    assert all(i == j for i, j in zip(vdj.data_names, vdj.data.index))
    assert all(i == j for i, j in zip(vdj.metadata_names, vdj.metadata.index))
