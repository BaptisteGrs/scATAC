#!/usr/bin/bash

for file in $1*;do
    file2=${file##*/}
    name=${file2%.bam}
    qsub -N "${name}_INDEX" Projects/DNAme/Scripts/index_bam_cluster.sh $file
done