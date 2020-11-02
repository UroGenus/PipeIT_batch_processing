#!/usr/bin/python3

import pandas as pd
import os
from os import path
import argparse
import sys

def main():
	parser = argparse.ArgumentParser( description="Collect PipeIT results in one file" )
	parser.add_argument( "-i", help="Input folder with PipeIT results", default = 'PipeIT/results')
	parser.add_argument( "-x", help="Ion Xpress Barcodes xlsx file with columns 'Sample Name', 'Normalize by', and 'Ion Xpress Barcode'", default = '/storage/research/dbmr_urology/Prostate_PDO/20200716_prostate_panel_sequencing.xlsx')
	parser.add_argument( "-o", help="Output .tsv file, default stdout", default = sys.stdout)

	pa = parser.parse_args()

	barcodes_df = pd.read_excel(pa.x, index_col=0, skiprows=[0])

	results_df = pd.DataFrame()
	for d in os.listdir(pa.i):
		if barcodes_df[barcodes_df['Sample Name'] == d].shape[0] == 0: continue
	
		results_file = path.join(pa.i, '%s/%s.PipeIT.tsv') % (d, d)
		if not path.isfile(results_file):
			print('No file', results_file)
			continue
		r_df = pd.read_csv(results_file, delimiter = '\t')
		if r_df.shape[0] == 0: continue
		norm_name = barcodes_df[barcodes_df['Sample Name'] == d]['Normalize by'].values[0]
		r_df.insert(loc=0, column='Sample Name', value=[d] * r_df.shape[0])
		r_df.insert(loc=1, column='Normalized by', value=[norm_name] * r_df.shape[0])
		results_df = results_df.append(r_df)

	results_df.sort_values(by = 'Sample Name', inplace = True)
	results_df.to_csv(pa.o, sep = '\t', index = False)

if __name__ == "__main__":
    main()
