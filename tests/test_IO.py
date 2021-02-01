#!/usr/bin/env python
# basic requirements for test data
import sys
import dandelion as ddl
import pandas as pd
import requests
from io import StringIO

def test_init():
	file = 'https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_B_1k_multi_5gex_b/sc5p_v2_hs_B_1k_multi_5gex_b_vdj_b_airr_rearrangement.tsv'
	r = requests.get(file)
	test_data = pd.read_csv(StringIO(r.text), sep = '\t')
	test_data['locus'] = ['IGH' if 'IGH' in i else 'IGK' if 'IGK' in i else 'IGL' if 'IGL' in i else None for i in test_data.v_call]
	test_data['sample_id'] = '10X'
	test_ddl = ddl.Dandelion(test_data)	
	test_ddl.write_h5('test.h5', compression = 'bzip2')
	test_ddl.write_pkl('test.pkl.pbz2')
	test2 = ddl.read_h5('test.h5')
	test3 = ddl.read_h5('test.pkl.pbz2')
	print(test2)
	print(test3)

if __name__ == '__main__':
	test_init()
