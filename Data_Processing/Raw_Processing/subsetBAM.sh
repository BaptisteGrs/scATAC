#!/usr/bin/bash
#$ -cwd
#$ -V
#$ -l h_rss=10G
#$ -pe smp 2

module load samtools

bam_file=$1
header=$2
filter=$3
path=$4

sample=${bam_file%.possorted_bam.bam}
cluster=${filter%_Sample_*}
export_file_bodysam="${path}${sample##*/}_${cluster##*/}_SAM_body"
export_file_sam="${path}${sample##*/}_${cluster##*/}.sam"
export_file="${path}${sample##*/}_${cluster##*/}.bam"
# Filter alignments 
samtools view $bam_file | LC_ALL=C grep -F -f $3 > $export_file_bodysam

# Combine header and body 
cat $2 $export_file_bodysam > $export_file_sam 

# Convert filtered.sam to BAm format
samtools view -b  $export_file_sam > $export_file