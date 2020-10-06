#!/usr/bin/python3

import pandas as pd
import os
from os import path
import argparse

def gen_str_code(code):
	if code < 10: return '00'+str(code)

	if code < 100: return '0'+str(code)

	return str(code)

def main():
	parser = argparse.ArgumentParser( description="Run PipeIT batch processing on Ubelix" )
	parser.add_argument( "-b", help="Folder with .bam data", default = '/storage/research/dbmr_urology/Prostate_PDO')
	parser.add_argument( "-i", help="PipeIT image file", default = '/home/ubelix/dbmr/ko20g613/PipeIT_1.2.13.img')
	parser.add_argument( "-t", help="Input bam folder", default = '/storage/research/dbmr_urology/Prostate_PDO/bam/')
	parser.add_argument( "-e", help="Target panel bed file", default = '/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed')
	parser.add_argument( "-x", help="Ion Xpress Barcodes xlsx file with columns columns 'Sample Name', 'Normalize by', and 'Ion Xpress Barcode'", default = '/storage/research/dbmr_urology/Prostate_PDO/20200716_prostate_panel_sequencing.xlsx')
	parser.add_argument( "-m", help="Email to report when jobs are done", required = False)
	parser.add_argument( "-s", help="snpEff jar file location", default = '/home/ubelix/dbmr/ko20g613/snpEff/SnpSift.jar')

	pa = parser.parse_args()

	# GENERATE ARGUMENTS FOR THE PIPELINE
	barcodes_df = pd.read_excel(pa.x, index_col=0, skiprows=[0])
	n_jobs = 0
	with open ('args.txt', 'w') as args_file:
		for norm, rows in barcodes_df.groupby('Normalize by'):
			norm_code = barcodes_df[barcodes_df['Sample Name'] == norm]['Ion Xpress Barcode'].values[0]
			norm_file = path.join(pa.t, "IonXpress_%s.bam") % (gen_str_code(norm_code))
			for i, row in rows.iterrows():
				n_jobs += 1
				tumour_code = row['Ion Xpress Barcode']
				tumour_file = path.join(pa.t, "IonXpress_%s.bam") % (gen_str_code(tumour_code))
				output_dir = path.join('PipeIT/results/', row['Sample Name'])
				vcf_file = '%s/%s.PipeIT.vcf'	% (output_dir, row['Sample Name'])	
				output_file = '%s/%s.PipeIT.tsv' % (output_dir, row['Sample Name'])

				args_file.write('%s\t%s\t%s\t%s\t%s\n' % (tumour_file, norm_file, row['Sample Name'], vcf_file, output_file))


	# CREATE A BASH EXECUTABLE FILE
	with open ('jobs.sh', 'w') as jsh:
		jsh.write('''\
#!/bin/bash

#SBATCH --mail-type=end,fail
#SBATCH --partition=all
#SBATCH --time=00:50:00    # Each task takes max 50 minutes
#SBATCH --mem-per-cpu=4G   # Each task uses max 4G of memory
#SBATCH --output=log.txt
#SBATCH --error=err.txt
''')
		jsh.write('#SBATCH --array=1-%i%%20  # Submit %i tasks with task ID 1,2,...,%i. Run max 20 tasks concurrently\n' % (n_jobs, n_jobs, n_jobs))
		if pa.m:
			jsh.write('\n#SBATCH --mail-user=%s\n') % (pa.m)

		jsh.write('\nmodule load Java/11.0.2\n\n')
		jsh.write('''\
param_store=args.txt

tumour_file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')  
norm_file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')
sample_name=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3}')  
vcf_file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $4}')
output_file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $5}')
''')

		jsh.write('\nsingularity run -B %s %s -t $tumour_file -n $norm_file -e %s -o $sample_name && java -jar %s extractFields -s "," $vcf_file CHROM POS REF ALT ANN[*].GENE ANN[*].HGVS_P AF > $output_file\n' % (pa.b, pa.i, pa.e, pa.s))
		
		jsh.write('\nexit')
	os.chmod('jobs.sh', 0o777)

if __name__ == "__main__":
    main()

