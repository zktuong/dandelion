#!/usr/bin/env python
# basic requirements for test data
import sys
import os
import dandelion as ddl
import scanpy as sc
import pandas as pd
import requests
from io import StringIO
from numba.core.errors import NumbaWarning, NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import pytest
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

def test_IO():
    file = "https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_B_1k_multi_5gex_b/sc5p_v2_hs_B_1k_multi_5gex_b_vdj_b_airr_rearrangement.tsv"
    r = requests.get(file)
    test_data = pd.read_csv(StringIO(r.text), sep="\t")
    test_data["locus"] = [
        "IGH" if "IGH" in i else "IGK" if "IGK" in i else "IGL" if "IGL" in i else None
        for i in test_data.v_call
    ]
    test_data["umi_count"] = test_data["duplicate_count"]
    test_data["sample_id"] = "test"
    test_ddl = ddl.Dandelion(test_data)
    test_ddl.write_h5("test/test.h5", compression="bzip2")
    test_ddl.write_pkl("test/test.pkl.pbz2")
    test = ddl.read_h5("test/test.h5")
    _ = ddl.read_pkl("test/test.pkl.pbz2")
    print(test)

@pytest.mark.filterwarnings('ignore::NumbaWarning')
def test_scanpy():
    scfile = "https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_B_1k_multi_5gex_b/sc5p_v2_hs_B_1k_multi_5gex_b_count_filtered_feature_bc_matrix.h5"
    r = requests.get(scfile)
    open("test/sctest.h5", "wb").write(r.content)
    adata = sc.read_10x_h5("test/sctest.h5")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata)
    adata.write("test/sctest.h5ad", compression="gzip")
    print(adata)


def test_filter():
    adata = sc.read_10x_h5("test/sctest.h5")
    test = ddl.read_h5("test/test.h5")
    adata.obs["filter_rna"] = False
    test, adata = ddl.pp.filter_bcr(test, adata)
    adata.write("test/sctest.h5ad", compression="gzip")
    test.write_h5("test/test.h5", compression="bzip2")


def test_update_metadata():
    test = ddl.read_h5("test/test.h5")
    ddl.update_metadata(test, "sequence_id")


def test_find_clones():
    test = ddl.read_h5("test/test.h5")
    ddl.tl.find_clones(test)
    test.write_h5("test/test.h5", compression="bzip2")


def test_generate_network():
    test = ddl.read_h5("test/test.h5")
    ddl.tl.generate_network(test, key="sequence_alignment")
    test.write_h5("test/test.h5", compression="bzip2")


def test_downsampling():
    test = ddl.read_h5("test/test.h5")
    test_downsample = ddl.tl.generate_network(
        test, key="sequence_alignment", downsample=100
    )
    print(test_downsample)


def test_transfer():
    test = ddl.read_h5("test/test.h5")
    adata = sc.read_h5ad("test/sctest.h5ad")
    ddl.tl.transfer(adata, test)
    adata.write("sctest.h5ad", compression="gzip")


def test_create_germlines():
    test = ddl.read_h5("test/test.h5")
    test.update_germline(germline="database/germlines/imgt/human/vdj/")
    ddl.pp.create_germlines(
        test,
        germline="database/germlines/imgt/human/vdj/",
        v_field="v_call",
        germ_types="dmask",
    )
    test.write_h5("test/test.h5", compression="bzip2")


def test_define_clones():
    test = ddl.read_h5("test/test.h5")
    ddl.pp.calculate_threshold(test, plot=False)
    ddl.tl.define_clones(test, key_added="changeo_clone_id")


def test_quantify_mutations():
    test = ddl.read_h5("test/test.h5")
    ddl.pp.quantify_mutations(test, germline_column="germline_alignment")


if __name__ == "__main__":
    test_airr()
    test_scanpy()
    test_filter()
    test_update_metadata()
    test_find_clones()
    test_generate_network()
    test_downsampling()
    test_transfer()
    test_create_germlines()
    test_define_clones()
    test_quantify_mutations()
