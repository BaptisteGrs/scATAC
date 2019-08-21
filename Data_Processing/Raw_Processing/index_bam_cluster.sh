#!/usr/bin/bash
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -l h_rss=8G
#$ -l h_vmem=8G
#$ -m ea
#$ -w e
#$ -M bgross@nygenome.org

module load samtools

samtools sort $1 "$1.sorted"
samtools index "$1.sorted.bam"