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
