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
	parser.add_argument( "-i", help="PipeIT image file", default = '/storage/homefs/ko20g613/PipeIT_2.0.0.img')
	parser.add_argument( "-t", help="Input bam folder", default = '/storage/research/dbmr_urology/Prostate_PDO/bam/')
	parser.add_argument( "-e", help="Target panel bed file", default = '/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed')
	parser.add_argument( "-x", help="Ion Xpress Barcodes xlsx file with columns columns 'Sample Name', 'Normalize by', and 'Ion Xpress Barcode'", default = '/storage/research/dbmr_urology/Prostate_PDO/20200716_prostate_panel_sequencing.xlsx')
	parser.add_argument( "-s", help="snpEff jar file location", default = '/storage/homefs/ko20g613/snpEff/SnpSift.jar')
	parser.add_argument( "-c", help="Annovar's database files folder", default = '/storage/research/dbmr_urology/Prostate_PDO/humandb')
	parser.add_argument( "-d", help="VCF file with the mutations found in a pool of normal samples", default = '/storage/research/dbmr_urology/Prostate_PDO/pon.tvc.vcf')
	parser.add_argument( "-m", help="Email to report when jobs are done (optional)", required = False)
	parser.add_argument( "-p", help="Separate character in the Ion Torrent .bam file names (default is '_')", default = '_')

	pa = parser.parse_args()

	# GENERATE ARGUMENTS FOR THE PIPELINE
	barcodes_df = pd.read_excel(pa.x, index_col=0, skiprows=[0])
	n_jobs_norm = 0
	n_jobs_nonorm = 0
	with open ('args_norm.txt', 'w') as args_norm_file:
		with open ('args_nonorm.txt', 'w') as args_nonorm_file:
			for i, row in barcodes_df[barcodes_df['Sample Type']!='Blood'].iterrows(): 
				tumour_code = row['Ion Xpress Barcode']
				tumour_file = path.join(pa.t, "IonXpress%s%s.bam") % (pa.p, gen_str_code(tumour_code))

				output_dir = path.join('PipeIT/results/', row['Sample Name'])
				vcf_file = '%s/%s.PipeIT.vcf'	% (output_dir, row['Sample Name'])	
				output_file = '%s/%s.PipeIT.tsv' % (output_dir, row['Sample Name'])
				
				norm = row['Normalize by']
				if pd.notna(norm):
					n_jobs_norm += 1
					norm_code = barcodes_df[barcodes_df['Sample Name'] == norm]['Ion Xpress Barcode'].values[0]
					norm_file = path.join(pa.t, "IonXpress%s%s.bam") % (pa.p, gen_str_code(norm_code))
					args_norm_file.write('%s\t%s\t%s\t%s\t%s\n' % (tumour_file, norm_file, row['Sample Name'], vcf_file, output_file))
				else:
					n_jobs_nonorm += 1					
					args_nonorm_file.write('%s\t%s\t%s\t%s\n' % (tumour_file, row['Sample Name'], vcf_file, output_file))


	# CREATE A BASH EXECUTABLE FILE FOR SAMPLES NORMALIZED WITH BLOOD
	with open ('jobs_norm.sh', 'w') as jsh:
		jsh.write('''\
#!/bin/bash

#SBATCH --mail-type=end,fail
#SBATCH --partition=epyc2
#SBATCH --time=04:00:00    # Each task takes max 4 hours
#SBATCH --mem-per-cpu=4G   # Each task uses max 4G of memory
#SBATCH --output=log_norm.txt
#SBATCH --error=err_norm.txt
''')
		jsh.write('#SBATCH --array=1-%i%%20  # Submit %i tasks with task ID 1,2,...,%i. Run max 20 tasks concurrently\n' % (n_jobs_norm, n_jobs_norm, n_jobs_norm))
		if pa.m:
			jsh.write('\n#SBATCH --mail-user=%s\n') % (pa.m)

		jsh.write('\nmodule load Java/11.0.2\n\n')
		jsh.write('''\
param_store=args_norm.txt

tumour_file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')  
norm_file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')
sample_name=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3}')  
vcf_file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $4}')
output_file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $5}')
''')

		jsh.write('\nsingularity run -B %s %s -t $tumour_file -n $norm_file -e %s -o $sample_name && java -jar %s extractFields -s "," $vcf_file CHROM POS REF ALT ANN[*].GENE ANN[*].GENEID ANN[*].FEATUREID ANN[*].HGVS_P AF > $output_file\n' % (pa.b, pa.i, pa.e, pa.s))
		
		jsh.write('\nexit')


	# CREATE A BASH EXECUTABLE FILE FOR SAMPLES NORMALIZED WITH SAMPLE POOL
	with open ('jobs_nonorm.sh', 'w') as jsh:
		jsh.write('''\
#!/bin/bash

#SBATCH --mail-type=end,fail
#SBATCH --partition=epyc2
#SBATCH --time=06:00:00    # Each task takes max 6 hours
#SBATCH --mem-per-cpu=4G   # Each task uses max 4G of memory
#SBATCH --output=log_nonorm.txt
#SBATCH --error=err_nonorm.txt
''')
		jsh.write('#SBATCH --array=1-%i%%20  # Submit %i tasks with task ID 1,2,...,%i. Run max 20 tasks concurrently\n' % (n_jobs_nonorm, n_jobs_nonorm, n_jobs_nonorm))
		if pa.m:
			jsh.write('\n#SBATCH --mail-user=%s\n') % (pa.m)

		jsh.write('\nmodule load Java/11.0.2\n\n')
		jsh.write('''\
param_store=args_nonorm.txt

tumour_file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')  
sample_name=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')  
vcf_file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3}')
output_file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $4}')
''')

		# PipeIT_1.2.13 version
		#jsh.write('\nsingularity run -B %s %s -t $tumour_file -e %s -c %s -d %s -r 4 -o $sample_name && java -jar %s extractFields -s "," $vcf_file CHROM POS REF ALT ANN[*].GENE ANN[*].GENEID ANN[*].FEATUREID ANN[*].HGVS_P AF > $output_file\n' % (pa.b, pa.i, pa.e, pa.c, pa.d, pa.s))
		
		# PipeIT_2.0.0 version
		jsh.write('\nsingularity run -B %s %s -t $tumour_file -e %s -c %s -o $sample_name && java -jar %s extractFields -s "," $vcf_file CHROM POS REF ALT ANN[*].GENE ANN[*].GENEID ANN[*].FEATUREID ANN[*].HGVS_P AF > $output_file\n' % (pa.b, pa.i, pa.e, pa.c, pa.d, pa.s))

		jsh.write('\nexit')

	os.chmod('jobs_norm.sh', 0o777)
	os.chmod('jobs_nonorm.sh', 0o777)

if __name__ == "__main__":
    main()

