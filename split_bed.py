#!/usr/bin/python3

import pandas as pd
import os
from os import path
import argparse
import sys
import numpy as np

def main():
	parser = argparse.ArgumentParser( description="Split bed file" )
	parser.add_argument( "-i", help="Bed file", default = '')

	pa = parser.parse_args()

	bed_df = pd.read_csv(pa.i, delimiter = '\t', header=None, usecols = [0, 1, 2], names = ['chr', 'start', 'end'])

	n = int(np.ceil(bed_df.shape[0] / 1000))

	for i in range(0, n):
		start = i*1000
		end = min(bed_df.shape[0], (i+1)*1000)
		rows = bed_df.iloc[start:end]
		rows.to_csv(pa.i+'_'+str(i)+'.bed', sep = '\t', index = False, header = False)
		print(rows.shape[0], start, end)

if __name__ == "__main__":
    main()
