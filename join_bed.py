#!/usr/bin/python3

import pandas as pd
import os
from os import path
import argparse
import sys

def main():
	parser = argparse.ArgumentParser( description="Join bed files" )
	parser.add_argument( "-i1", help="Original bed file", default = '')
	parser.add_argument( "-i2", help="Bed file with gene names", default = '')
	parser.add_argument( "-o1", help="Output region bed file", default = '')
	parser.add_argument( "-o2", help="Output gene bed file", default = '')

	pa = parser.parse_args()

	region_df = pd.read_csv(pa.i1, delimiter = '\t', header=None, usecols = [0, 1, 2, 3, 5], names = ['chr', 'region_start', 'region_end', 'region_id', 'gene_pool'])
	gene_df = pd.read_csv(pa.i2, delimiter = '\t', header=None, usecols = [0, 1, 2, 3, 4, 5], names = ['chr', 'chr_start', 'chr_end', 'geneName', 'knownToEnsembl.name', 'knownToEnsembl.value'])
	gene_df = gene_df.drop_duplicates()

	print('Number of gene names =', len(gene_df['geneName'].unique()))
	print('Number of lines in the genes file =',gene_df.shape[0])
	gene_df.to_csv(pa.o2, sep = '\t', index = False)
	
	for i, row in region_df.iterrows():
		gene_rows = gene_df[(gene_df.chr_start <= row.region_start) & (gene_df.chr_end >= row.region_end)]
		region_df.loc[i, 'geneNames'] = ','.join(gene_rows['geneName'].fillna('').values)
		region_df.loc[i, 'knownToEnsembl.name'] = ','.join(gene_rows['knownToEnsembl.name'].fillna('').values)
		region_df.loc[i, 'knownToEnsembl.value'] = ','.join(gene_rows['knownToEnsembl.value'].fillna('').values)
	
	region_df.to_csv(pa.o1, sep = '\t', index = False)

if __name__ == "__main__":
    main()
