# scATAC pipeline : clustering and motif analysis

Pipeline to analyze scATAC-seq data using cisTopic, HOMER and chromVAR

## A few introduction notes 

The data has been produced thanks to 10X Genomics scATAC-seq technology. 
It consists of 6 different samples from mice's bone marrow. 3 genotypes were studied : control, DMT3A KO and TET2 KO. 

Sequence files for each of the samples were aligned and processed thanks to `cellranger-atac`. The results were then aggregated using `cellranger aggr`. Identifiers were added to the end of the barcodes to distinguish between samples (-1 and -2 for DNMT3A KO samples, -3 and -4 for TET2 KO samples and -5 -6 for WT samples)

cellRanger provides BAM files for each of the samples, a merged matrix of peaks * cells and a fragments.tsv.gz file.
You can learn more about the outputs of cellRanger here (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices)

I used cisTopic for dimensionality reduction (https://github.com/aertslab/cisTopic/) but other methods exist and perform well. For future analysis, I would recommand also taking a look at : 
- https://github.com/r3fang/SnapATAC 
- https://github.com/timoast/signac

They seem to produce accurate clustering for reasonable running times (cisTopic may be a bit longer). 

I'd also **highly** recommand taking a look at these 2 articles that benchmark scATAC analysis methods : 
- https://www.biorxiv.org/content/10.1101/739011v1
- http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/

Overall, cisTopic is still listed as one of the best ways to analyze scATAC and also provide an interesting integration of the  AUCell library(https://github.com/aertslab/AUCell), that allows to look at the enrichment of epigenetic signatures in the cells. 

## Repository organisation 

