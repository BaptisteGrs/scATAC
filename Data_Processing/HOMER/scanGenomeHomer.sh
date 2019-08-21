#!/bin/bash
#$ -V
path='Projects/DNAme/data/HOCOMOCO_mouse/'

for file in $path*;do
	scanMotifGenomeWide.pl $file mm10 -bed > ${file%.motif}'.sites.bed' 
done
 
