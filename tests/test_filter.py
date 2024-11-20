#!/usr/bin/env python
import json
import os
import dandelion as ddl
import pytest


@pytest.mark.usefixtures("airr_generic")
def test_generic_check(airr_generic):
    """test data loading and filtering"""
    tmp = ddl.Dandelion(airr_generic)
    assert tmp.metadata.shape[0] == 45
    assert tmp.data.shape[0] == airr_generic.shape[0]

    tmp2 = ddl.pp.check_contigs(tmp, productive_only=False)
    assert tmp2.metadata.shape[0] == 44
    assert (
        tmp2.data.shape[0] != tmp.data.shape[0]
    )  # this is because filter_extra is True by default

    tmp2 = ddl.pp.check_contigs(tmp, productive_only=False, library_type="ig")
    assert tmp2.metadata.shape[0] == 25
    assert tmp2.data.shape[0] != tmp.data.shape[0]
    assert tmp2.data.shape[0] == 65

    tmp2 = ddl.pp.check_contigs(tmp)
    assert tmp2.metadata.shape[0] == 43
    assert tmp2.data.shape[0] != tmp.data.shape[0]

    ddl.tl.find_clones(tmp2, identity={"tr-ab": 1})
    assert "clone_id" in tmp2.data
    assert "clone_id" in tmp2.metadata
    assert not tmp2.metadata.clone_id.empty

    ddl.tl.generate_network(tmp2, key="junction_aa", compute_layout=False)
    assert tmp2.graph is not None


@pytest.mark.usefixtures("airr_generic")
def test_check_keep_extra(airr_generic):
    """test keep extra option"""
    tmp = ddl.Dandelion(airr_generic)
    assert tmp.metadata.shape[0] == 45
    assert tmp.data.shape[0] == airr_generic.shape[0]
    tmp2 = ddl.pp.check_contigs(tmp)
    tmp3 = ddl.pp.check_contigs(tmp, filter_extra=False)
    assert tmp3.metadata.shape[0] == 43
    assert tmp2.data.shape[0] != tmp.data.shape[0]
    assert tmp3.data.shape[0] == tmp.data.shape[0]
    assert tmp3.data.extra.value_counts()["T"] == 14


@pytest.mark.usefixtures("airr_generic")
def test_check_remove_ambiguous(airr_generic):
    """test remove ambiguous option"""
    tmp = ddl.Dandelion(airr_generic)
    assert tmp.metadata.shape[0] == 45
    assert tmp.data.shape[0] == airr_generic.shape[0]
    tmp2 = ddl.pp.check_contigs(tmp)
    tmp3 = ddl.pp.check_contigs(tmp, filter_ambiguous=True)
    assert tmp3.metadata.shape[0] == 43
    assert tmp2.data.shape[0] != tmp.data.shape[0]
    assert tmp3.data.shape[0] != tmp2.data.shape[0]
    assert tmp2.data.ambiguous.value_counts()["T"] == 37
    assert tmp2.data.ambiguous.value_counts()["F"] == 79
    assert all(tmp3.data.ambiguous == "F")
