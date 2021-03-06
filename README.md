# PipeIT batch processing

This repository contains scripts for batch processing of raw .bam files using [PipeIT](https://github.com/ckynlab/PipeIT) pipeline. [PipeIT](https://github.com/ckynlab/PipeIT) is a somatic variant calling pipeline specific for Ion Torrent sequencing data.

## Batch processing on UBELIX

### Input data structure
- Folder with .bam files, so that each file has a name `IonXpress_<code>.bam`, see examples in `/storage/research/dbmr_urology/Prostate_PDO/bam`
- xlsx file describing the .bam files with the columns 'Sample Name', 'Normalize by', and 'Ion Xpress Barcode' with the same structure as the file `/storage/research/dbmr_urology/Prostate_PDO/20200716_prostate_panel_sequencing.xlsx`
- .bed file with the target gene panel, see example in `/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed`

### Software to be downloaded
- `PipeIT_<version>.img` is PipeIT singularity image that can be downloaded from PipeIT laboratory's website: [http://oncogenomicslab.org/software-downloads](http://oncogenomicslab.org/software-downloads)
- SnpEff & SnpSift is a genomic variant annotations and functional effect prediction toolbox that can be installed as described at [https://pcingola.github.io/SnpEff/download/](https://pcingola.github.io/SnpEff/download/)

### Processing
Run the following command with the parameters:
```
python3 PipeIT_batch_run_ubelix.py
  -b B        Folder with .bam files
  -i I        PipeIT image file
  -t T        Input bam folder
  -e E        Target panel bed file
  -x X        Ion Xpress Barcodes xlsx file with columns columns 'Sample Name', 'Normalize by', and 'Ion Xpress Barcode'
  -s S        snpEff jar file
  -c C        Annovar's database files folder
  -d D        VCF file for the mutations found in a pool of normal samples
  -m M        Email to report when jobs are done (optional)
  -p P        Separate character in the Ion Torrent .bam file names (default is '_')
```
For example
```
python3 PipeIT_batch_run_ubelix.py -i /home/ubelix/dbmr/ko20g613/PipeIT_1.2.13.img -b /storage/research/dbmr_urology/Prostate_PDO /home/ubelix/dbmr/ko20g613/PipeIT_1.2.13.img -t /storage/research/dbmr_urology/Prostate_PDO/bam -e /storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed -x /storage/research/dbmr_urology/Prostate_PDO/20200716_prostate_panel_sequencing.xlsx -s /home/ubelix/dbmr/ko20g613/snpEff/SnpSift.jar -c /storage/research/dbmr_urology/humandb -d /storage/research/dbmr_urology/Prostate_PDO/pon.tvc.vcf
```
This script will generate the `jobs_norm.sh` and `jobs_nonorm.sh` files in the same directory, that should be run as follows
```
sbatch jobs_norm.sh
sbatch jobs_nonorm.sh
```

Processing can be monitored by the command `squeue --user <user name>` as described [here](https://hpc-unibe-ch.github.io/user-guide/job-management/monitoring-jobs.html). Logs can be found in the files `log_(no)norm.txt` and `err_(no)norm.txt` in the same directory.

After the processing has finished, the output .vcf and .tsv files for each sample will be located in the folder `PipeIT/results/<sample name>`.
You can then collect all .tsv results in one file by running the following command
```
python3 collect_PipeIT_results.py
  -i I        Input folder with PipeIT results
  -x X        Ion Xpress Barcodes xlsx file with columns 'Sample Name', 'Normalize by', and 'Ion Xpress Barcode'
  -o O        Output .tsv file, default stdout

```
For example
```
python3 collect_PipeIT_results.py -i PipeIT/results -x /storage/research/dbmr_urology/Prostate_PDO/20200716_prostate_panel_sequencing.xlsx -o PDO_mutation_results.tsv
```

## Process one file

### With blood sample for normalization
To process one file with blood sample for normalization, run the following command
```
singularity run -B <path/to/folders/that/need/to/be/mounted> <path/to/PipeIT_<version>.img> -t path/to/tumor.bam -n path/to/normal.bam -e path/to/region.bed  -o <output_directory>
```
For example, for prostate PDO raw data run on UBELIX cluster
```
singularity run -B /storage/research/dbmr_urology/Prostate_PDO PipeIT_1.2.13.img -t /storage/research/dbmr_urology/Prostate_PDO/bam/IonXpress_003.bam -n /storage/research/dbmr_urology/Prostate_PDO/bam/IonXpress_010.bam -e /storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed -o test_3
```
where
- `/storage/research/dbmr_urology/Prostate_PDO` is the folder with prostate PDO raw data
- `PipeIT_1.2.13.img` is PipeIT singularity image that can be downloaded from PipeIT laboratory's website: [http://oncogenomicslab.org/software-downloads](http://oncogenomicslab.org/software-downloads)
- `/storage/research/dbmr_urology/Prostate_PDO/bam/IonXpress_003.bam` is the path to the tumour raw .bam file
- `/storage/research/dbmr_urology/Prostate_PDO/bam/IonXpress_001.bam`is the path to the blood raw .bam file used for normalizing
- `/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed` is the path to the .bed file with the targeted gene regions
- `test_3` is the name of the output folder that will be generated by PipeIT

### With normal sample pool for normalization
To process one file with normal sample pool for normalization, run the following command
```
singularity run -B <path/to/folders/that/need/to/be/mounted> <path/to/PipeIT_<version>.img> -t path/to/tumor.bam -e path/to/region.bed -c path/to/annovar/folder -d /path/to/mutations.vcf -o <output_directory>
```
For example, for prostate PDO raw data run on UBELIX cluster
```
singularity run -B /storage/research/dbmr_urology/Prostate_PDO /home/ubelix/dbmr/ko20g613/PipeIT_1.2.13.img -t /storage/research/dbmr_urology/Prostate_PDO/bam/IonXpress_024.bam  -e /storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed -c /storage/research/dbmr_urology/Prostate_PDO/humandb -d /storage/research/dbmr_urology/Prostate_PDO/pon.tvc.vcf  -r 4  -o test_24
```
where
- `/storage/research/dbmr_urology/Prostate_PDO/humandb` is the path to the Annovar's database folder that can be dowloaded with PipeIT as described [here](https://github.com/ckynlab/PipeIT)
- `/storage/research/dbmr_urology/Prostate_PDO/pon.tvc.vcf` is the path to the file with the mutations found in a pool of normal samples (can be generated with [this pipeline](https://github.com/charlottekyng/usb-modules-v2))

For more details and options, see [the PipeIT repository docu](https://github.com/ckynlab/PipeIT).

### Generate final output as .tsv

PipeIT generates several output files in the folder `PipeIT/results/<output_directory>`. The .vcf can be found under `PipeIT/results/<output_directory>/<output_directory>.PipeIT.vcf`. To extract mutation information from this file in .tsv format, install [SnpEff & SnpSift](https://pcingola.github.io/SnpEff/download/) and run the following commands

On UBELIX, first the java module has to be loaded as follows: `module load <java module name>`. For example, `module load Java/11.0.2`. Then run
```
java -jar path/to/SnpSift.jar extractFields -s "," PipeIT/results/<output_directory>/<output_directory>.PipeIT.vcf CHROM POS REF ALT ANN[*].GENE ANN[*].GENEID ANN[*].FEATUREID ANN[*].HGVS_P AF > output.tsv
```

To extract more information from .vcf, check out [documentation of SnpSift](https://pcingola.github.io/SnpEff/ss_extractfields/#example-1-extracting-chromosome-position-id-and-allele-frequency).

## Retrieve Gene Symbol IDs from the .bed file

If you want to associate the official gene symbol to the regions overlapping the human genome indicated in the .bed file, you can 
navigate to [Genome Browser Gateway](https://genome.ucsc.edu/cgi-bin/hgGateway) and select table browser from the tools menu.

Complete the fields in the table browser with [these options](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=915189327_SVlXMVfDA3Fea7LjM0AaKepVBWlP&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=knownGene&hgta_table=knownCanonical&hgta_regionType=userRegions&position=chrX%3A15%2C578%2C261-15%2C621%2C068&hgta_outputType=selectedFields&hgta_outFileName=output.test.01), select "define region" and click "change. Upload the .bed file corresponding to the IonTorrnet panel and click "submit". **Note:** the input file should be <= 1000 lines (split the original .bed file using [this script](https://github.com/UroGenus/PipeIT_batch_processing/blob/main/split_bed.py)). 

Select the following fields from the get output view:

- **Select Fields from hg19.knownCanonical** : chrom, chromstart, chromend
- **hg19.kgXref fields** : genSymbol
- **hg19.knownToEnsembl fields** : primary ID, associated ID, label, quantity etc.
- **Linked Tables** : hg19 knownToEnsembl

Go back to the top of the page and click "get output" again, before saving the file.

The file obtained will have 6 columns:
- chromosome name
- chromosome start
- chromosome end
- gene symbol
- knownToEnsembl.name
- knownToEnsembl.value

Gene information can be then added to the original .bed file using [this script](https://github.com/UroGenus/PipeIT_batch_processing/blob/main/join_bed.py).
