#!/usr/bin/bash
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -l h_rss=20G
#$ -l h_vmem=20G
#$ -m ea
#$ -w e
#$ -M bgross@nygenome.org

module load samtools

var=$(ls $2* | grep -E "$3.*$4|$4.*$3")
echo $var
samtools merge -f $1 $var