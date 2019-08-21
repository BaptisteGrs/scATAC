#!/usr/bin/bash

path=$1
files=$(ls $1 | grep "[^bai]$" | grep .sorted.bam)
echo $files
for file in $files;do
    name=${file%.bam.sorted.bam}
    echo $file
    echo $name
    qsub -N "${name}_Normalize" Projects/DNAme/Scripts/bam_normalize.sh $path$file $path$name".bw"
done

