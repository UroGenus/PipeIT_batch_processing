module load vital-it/7
module load UHTS/Analysis/samtools/1.10

for BAM in $1/*.bam
do
	echo $BAM
	filename=$(basename -- "$BAM")
	filename="${filename%.*}"
	samtools view -H $BAM | sed "s/SM:.*$/SM:$filename/" | samtools reheader - $BAM > tmp.bam
	mv tmp.bam $BAM
done
