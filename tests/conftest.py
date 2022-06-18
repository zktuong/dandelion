#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2022-06-18 11:41:31
# @Last Modified by:   Kelvin
# @Last Modified time: 2022-06-18 15:17:00
"""conftest module."""
from .fixtures.fixtures import (
    create_testfolder,
    database_paths,
    processed_files,
    processed_files_tr,
    dummy_adata,
    dummy_adata2,
    dummy_adata_cr6,
    dummy_adata_tr,
    fasta_10x,
    annotation_10x,
    annotation_10x_cr6,
    fasta_10x_cr6,
    airr_10x,
    airr_reannotated,
    airr_reannotated2,
    fasta_10x_tr1,
    fasta_10x_tr2,
    annotation_10x_tr1,
    json_10x_cr6,
    annotation_10x_tr2,
    airr_travdv,
    fasta_10x_travdv,
    annotation_10x_travdv,
    dummy_adata_travdv,
)
from .fixtures.fixtures_mouse import (
    database_paths_mouse,
    dummy_adata_mouse,
    fasta_10x_mouse,
    annotation_10x_mouse,
    balbc_ighg_primers,
)
