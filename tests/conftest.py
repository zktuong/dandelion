#!/usr/bin/env python
# @Author: Kelvin
"""conftest module."""
from .fixtures.fixtures import (
    airr_10x,
    airr_generic,
    airr_reannotated,
    airr_reannotated2,
    airr_travdv,
    annotation_10x,
    annotation_10x_cr6,
    annotation_10x_tr1,
    annotation_10x_tr2,
    annotation_10x_travdv,
    create_testfolder,
    database_paths,
    dummy_adata,
    dummy_adata2,
    dummy_adata_cr6,
    dummy_adata_tr,
    dummy_adata_travdv,
    fasta_10x,
    fasta_10x_cr6,
    fasta_10x_tr1,
    fasta_10x_tr2,
    fasta_10x_travdv,
    json_10x_cr6,
    processed_files,
    processed_files_tr,
)
from .fixtures.fixtures_mouse import (
    annotation_10x_mouse,
    balbc_ighg_primers,
    database_paths_mouse,
    dummy_adata_mouse,
    fasta_10x_mouse,
)
