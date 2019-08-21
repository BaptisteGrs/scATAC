#!/bin/bash
#$ -V
#$ -pe smp 2
#$ -m ea
#$ -M bgross@nygenome.org

file=$1
out=$3
findMotifsGenome.pl $file mm10 $out/ -p 4 -size 200 -mask -bg $2 -mknown $4 -mcheck $4

