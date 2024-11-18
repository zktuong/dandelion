#!/usr/bin/env python
import dandelion as ddl
import pytest


@pytest.mark.usefixtures("airr_generic")
def test_query(airr_generic):
    """test query and update_metadata functions"""
    vdj = ddl.Dandelion(airr_generic)
    ddl.pp.check_contigs(vdj)
    vdj.update_metadata(retrieve="umi_count", retrieve_mode="split and sum")
    vdj.update_metadata(retrieve="umi_count", retrieve_mode="sum")
    vdj.update_metadata(retrieve="umi_count", retrieve_mode="average")
    vdj.update_metadata(retrieve="np2_length", retrieve_mode="split and sum")
    vdj.update_metadata(retrieve="np2_length", retrieve_mode="average")
    vdj.update_metadata(retrieve="np2_length", retrieve_mode="sum")

    vdj.update_metadata(
        retrieve="junction_aa",
        retrieve_mode="split and unique only",
        by_celltype=True,
    )
    vdj.update_metadata(
        retrieve="junction_aa",
        retrieve_mode="merge and unique only",
        by_celltype=True,
    )
    vdj.update_metadata(
        retrieve="junction_aa", retrieve_mode="merge", by_celltype=True
    )
    vdj.update_metadata(
        retrieve="junction_aa", retrieve_mode="split", by_celltype=True
    )
    vdj.update_metadata(
        retrieve="np2_length",
        retrieve_mode="split and average",
        by_celltype=True,
    )
    vdj.update_metadata(
        retrieve="np2_length", retrieve_mode="sum", by_celltype=True
    )
    vdj.update_metadata(
        retrieve="np2_length", retrieve_mode="average", by_celltype=True
    )
    vdj.update_metadata(
        retrieve="np2_length", retrieve_mode="split", by_celltype=True
    )
