#!/usr/bin/bash
#$-V
#$-cwd

module load samtools

path_to_XX=$1

bam_file=$(ls ${path_to_XX}/outs/ | grep possorted_bam.bam | grep -v .bai | grep -v .sam)
echo Bam file found : $bam_file

## First step is to modify the possorted bam file from 10X to have the cell barcode at the beginning of each read
echo Creating new bam file...
# Extract header
samtools view ${path_to_XX}/outs/${bam_file} -H > ${path_to_XX}/outs/${bam_file}.header.sam
samtools view ${path_to_XX}/outs/${bam_file} | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' > ${path_to_XX}/outs/temp.bam

cat ${path_to_XX}/outs/${bam_file}.header.sam ${path_to_XX}/outs/temp.bam | samtools view -bS - > ${path_to_XX}/outs/possorted.snap.bam

echo Deleting temporary file...
rm ${path_to_XX}/outs/temp.bam

# Sort the file by read name
echo Sorting new bam file by read name... 
samtools sort -n -@ 10 -m 1G ${path_to_XX}/outs/possorted.snap.bam ${path_to_XX}/outs/possorted.snap.nsrt

echo Creating snap object...
/gpfs/commons/home/bgross/anaconda3/envs/bioR/bin/snaptools snap-pre  \
                                                         --input-file=${path_to_XX}/outs/possorted.snap.nsrt.bam \
                                                         --output-snap=${path_to_XX}/outs/possorted.snap  \
                                                         --genome-name=mm10  \
                                                         --genome-size=/gpfs/commons/home/bgross/Projects/DNAme/data/mm10.chrom.sizes  \
                                                         --min-mapq=30  \
                                                         --min-flen=50  \
                                                         --max-flen=1000  \
                                                         --keep-chrm=TRUE  \
                                                         --keep-single=FALSE  \
                                                         --keep-secondary=False  \
                                                         --overwrite=True  \
                                                         --max-num=20000  \
                                                         --min-cov=500  \
                                                         --verbose=True

echo Done!                                                        
