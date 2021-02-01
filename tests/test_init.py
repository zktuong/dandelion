#!/usr/bin/env python
# basic requirements for test data
import sys
import os
import dandelion as ddl
import scanpy as sc
import pandas as pd
import requests
from io import StringIO

def test_IO():
	file = 'https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_B_1k_multi_5gex_b/sc5p_v2_hs_B_1k_multi_5gex_b_vdj_b_airr_rearrangement.tsv'
	r = requests.get(file)
	test_data = pd.read_csv(StringIO(r.text), sep = '\t')
	test_data['locus'] = ['IGH' if 'IGH' in i else 'IGK' if 'IGK' in i else 'IGL' if 'IGL' in i else None for i in test_data.v_call]
	test_data['umi_count'] = test_data['duiplicate_count']
	test_data['sample_id'] = '10X'
	test_ddl = ddl.Dandelion(test_data)
	test_ddl.write_h5('test.h5', compression = 'bzip2')
	test_ddl.write_pkl('test.pkl.pbz2')
	test = ddl.read_h5('test.h5')
	_ = ddl.read_pkl('test.pkl.pbz2')
	print(test)

def test_scanpy():
	scfile = 'https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_B_1k_multi_5gex_b/sc5p_v2_hs_B_1k_multi_5gex_b_count_filtered_feature_bc_matrix.h5' 
	r = requests.get(scfile)
	open('sctest.h5', 'wb').write(r.content)
	adata = sc.read_10x_h5('sctest.h5')
	print(adata)

def test_filter():
	adata = sc.read_10x_h5('sctest.h5')
	test = ddl.read_h5('test.h5')
	adata.obs['filter_rna'] = False
	test, adata = ddl.pp.filter_bcr(test, adata)
	adata.write('sctest.h5ad', compression = 'gzip')
	test.write_h5('test.h5', compression = 'bzip2')

def test_update_metadata():
	test = ddl.read_h5('test.h5')
	ddl.update_metadata(test, 'sequence_id')

def test_find_clones():
	test = ddl.read_h5('test.h5')
	ddl.tl.find_clones(test)

def test_generate_network():
	test = ddl.read_h5('test.h5')
	ddl.tl.generate_network(test, key = 'sequence_alignment')

def test_downsampling():
	test = ddl.read_h5('test.h5')
	ddl.tl.generate_network(test, key = 'sequence_alignment')
	test_downsample = ddl.tl.generate_network(test, key = 'sequence_alignment', downsample = 500)
	print(test_downsample)

if __name__ == '__main__':
	test_IO()
	test_scanpy()
	test_filter()
	test_update_metadata()
	test_find_clones()
	test_generate_network()
	test_downsampling()

