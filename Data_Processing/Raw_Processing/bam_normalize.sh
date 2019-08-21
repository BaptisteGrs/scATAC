#!/usr/bin/bash
#$ -cwd
#$ -V
#$ -l h_rss=100G
#$ -w e
#$ -p -1
#$ -M bgross@nygenome.org
#$ -m ea
#$ -pe smp 4

module load deeptools

bamCoverage --bam $1 --outFileName $2 \
            --normalizeUsingRPKM -bs 1 -p max -v 
        