#!/usr/bin/python3

import os
from os import path
import argparse
import sys

def main():
	parser = argparse.ArgumentParser( description="Collect PipeIT results in one file" )
	parser.add_argument( "-i", help="Input folder with PipeIT results", default = 'PipeIT/results')
	parser.add_argument( "-s", help="snpEff jar file location", default = '/home/ubelix/dbmr/ko20g613/snpEff/SnpSift.jar')

	pa = parser.parse_args()

	cmd = ('module load Java/11.0.2')
	os.system(cmd)
	for d in os.listdir(pa.i):
		results_file = path.join(pa.i, '%s/%s.PipeIT.vcf') % (d, d)
		if not path.isfile(results_file):
			print('No file', results_file)
			continue
		out_file = path.join(pa.i, '%s/%s.PipeIT.tsv') % (d, d)

		cmd = 'java -jar %s extractFields -s "," %s CHROM POS REF ALT ANN[*].GENE ANN[*].GENEID ANN[*].FEATUREID ANN[*].HGVS_P AF > %s' % (pa.s, results_file, out_file)
		os.system(cmd)
		

if __name__ == "__main__":
    main()
