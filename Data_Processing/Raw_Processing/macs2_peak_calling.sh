#!/usr/bin/bash
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -l h_rss=10G
#$ -l h_vmem=10G
#$ -w e

module load macs2

macs2 callpeak -t $1 -f BAM --broad --broad-cutoff 0.1 -m 2 100 --outdir $2 \
                -n $3