#!/bin/bash

clusters="Cluster_1 Cluster_2 Cluster_3 Cluster_4 Cluster_5 Cluster_6 Cluster_7"
genotypes="DNMT3A TET2 WT"

path=$1

for cluster in $clusters;do
    name="${path}Merged_${cluster}"
    qsub -N "Merged_${cluster}" Projects/DNAme/Scripts/samtools_merge.sh "${name}.bam" Projects/DNAme/analysis/KO_projection/ClusterBAM/Raw_files/Raw_BAM_files/ $cluster
    #for genotype in $genotypes;do
    #    name2="${path}Merged_${cluster}_${genotype}"
    #    qsub -N "Merged_${cluster}_${genotype}" Projects/DNAme/Scripts/samtools_merge.sh "${name2}.bam" Projects/DNAme/analysis/KO_projection/ClusterBAM/Bam_files/ $cluster $genotype
    #done
done