#!/usr/bin/env python
# basic requirements for test data
import sys
import os
from io import StringIO
import requests
import pandas as pd
import scanpy as sc
import dandelion as ddl


def test_setup():
    file = "https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_B_1k_multi_5gex_b/sc5p_v2_hs_B_1k_multi_5gex_b_vdj_b_airr_rearrangement.tsv"
    r = requests.get(file)
    test_data = pd.read_csv(StringIO(r.text), sep="\t")
    test_ddl = ddl.read_10x_airr(test_data)
    test_ddl.write_h5("tests/test.h5", compression="bzip2")
    test_ddl.write_pkl("tests/test.pkl.pbz2")
    test = ddl.read_h5("tests/test.h5")
    _ = ddl.read_pkl("tests/test.pkl.pbz2")
    scfile = "https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_B_1k_multi_5gex_b/sc5p_v2_hs_B_1k_multi_5gex_b_count_filtered_feature_bc_matrix.h5"
    r = requests.get(scfile)
    open("tests/sctest.h5", "wb").write(r.content)
    adata = sc.read_10x_h5("tests/sctest.h5")
    adata.write("tests/sctest.h5ad", compression="gzip")
    print(test)
    print(adata)


def test_filter():
    adata = sc.read_10x_h5("tests/sctest.h5")
    test = ddl.read_h5("tests/test.h5")
    adata.obs["filter_rna"] = False
    test, adata = ddl.pp.filter_bcr(test, adata)
    adata.write("tests/sctest.h5ad", compression="gzip")
    test.write_h5("tests/test.h5", compression="bzip2")
    print(test)


def test_update_metadata():
    test = ddl.read_h5("tests/test.h5")
    ddl.update_metadata(test, "sequence_id")
    print(test)


def test_find_clones():
    test = ddl.read_h5("tests/test.h5")
    ddl.tl.find_clones(test)
    test.write_h5("tests/test.h5", compression="bzip2")
    print(test)


def test_generate_network():
    test = ddl.read_h5("tests/test.h5")
    ddl.tl.generate_network(test, key="sequence_alignment")
    test.write_h5("tests/test.h5", compression="bzip2")
    print(test)


def test_downsampling():
    test = ddl.read_h5("tests/test.h5")
    test_downsample = ddl.tl.generate_network(
        test, key="sequence_alignment", downsample=100)
    print(test_downsample)


def test_transfer():
    test = ddl.read_h5("tests/test.h5")
    adata = sc.read_h5ad("tests/sctest.h5ad")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata)
    ddl.tl.transfer(adata, test)
    adata.write("tests/sctest.h5ad", compression="gzip")
    print(adata)


def test_create_germlines():
    test = ddl.read_h5("tests/test.h5")
    test.update_germline(germline="database/germlines/imgt/human/vdj/")
    ddl.pp.create_germlines(
        test, germline="database/germlines/imgt/human/vdj/", v_field="v_call", germ_types="dmask")
    test.write_h5("tests/test.h5", compression="bzip2")
    print(test)


def test_define_clones():
    test = ddl.read_h5("tests/test.h5")
    ddl.pp.calculate_threshold(test, plot=False)
    ddl.tl.define_clones(test, key_added="changeo_clone_id")
    print(test)


def test_quantify_mutations():
    test = ddl.read_h5("tests/test.h5")
    ddl.pp.quantify_mutations(test, germline_column="germline_alignment")
    print(test)


if __name__ == "__main__":
    test_setup()
    test_filter()
    test_update_metadata()
    test_find_clones()
    test_generate_network()
    test_downsampling()
    test_transfer()
    test_create_germlines()
    test_define_clones()
    test_quantify_mutations()
