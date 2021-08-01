#!/usr/bin/env python
import dandelion as ddl
import scanpy as sc


def test_recipe():
    adata = sc.datasets.pbmc3k()
    ddl.pp.external.recipe_scanpy_qc(adata)
    assert not adata.obs['filter_rna'].empty
    ddl.pp.external.recipe_scanpy_qc(adata, mito_cutoff = None)
    assert not adata.obs['gmm_pct_count_clusters_keep'].empty
