#!/usr/bin/env python
import pandas as pd
import sys

outname = sys.argv[1]
filenames = sys.argv[2]
option = sys.argv[3]

filenames = filenames.split(' ')
if option == '1':
	combined = pd.concat([pd.read_csv(f) for f in filenames])
	combined.to_csv(outname, index=None)
if option == '2':
	combined = pd.concat([pd.read_csv(f, sep = '\t') for f in filenames])
	combined.to_csv(outname, sep = '\t', index=None)
