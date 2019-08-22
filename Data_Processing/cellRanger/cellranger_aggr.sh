#!/bin/sh
#$-pe smp 10
#$-l h_rss=30G

module load cellranger-atac

cd /gpfs/commons/groups/landau_lab/Franco_Baptiste/

cellranger-atac aggr --id HematoProject_cellRanger \
                     --csv=aggregation.csv \
                     --normalize=none \
                     --reference=RawData/ref/refdata-cellranger-atac-mm10-1.1.0 \
                     --jobmode=local \
                     --localmem=30 \
                     --localcores=10 \
                     --nosecondary



