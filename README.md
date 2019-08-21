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

Overall, cisTopic is still listed as one of the best ways to analyze scATAC and also provide an interesting integration of the  AUCell library(https://github.com/aertslab/AUCell), that allows to look at the enrichment of epigenetic signatures in the cells (which is convenient for cluster annotation). 

## Repository organisation 

### `Data_Processing` folder

The `Data_Processing` folder contains scripts used to generate models from the scATAC data. 

#### cisTopic

The first step of the pipeline is to train an LDA model over the peaks-cells matrix from cellRanger.
The `train_cisTopic.R` script is using the `createcistopic` function from `utils/tools_factory.R` to initialize a cisTopic object, given specific filters (on barcodes and peaks coverage). 

Once the model is trained, clusters are created in the `Analysis/big_analysis.R` script by performing Louvain clustering on the cell-topic z-scores matrix (graph created using k=400 neighbors - k is to be chosen to have a consistent number of clusters and of suficient sizes).

Second step is to annotate the clusters. I've used `AUCell` to calculate a score for each epigenetic signatures of the cell types available in the ImmGen dataset (more about this in the `ImmGen_Preparation` paragraph). 
This is done by running the `signatures.R` script, that takes into argument a cisTopic model and the path to a directory containing a bed file for each signature you'd want to look at (in this case, epigenetic signatures by cell types). 

RCistarget can also be run to identify cistromes and for motif analysis, but it needs a feather file that was impossible for me to download at the time ('mm9-regions-9species.all_regions.mc9nr' on https://resources.aertslab.org/cistarget/) 

#### chromVAR

chromVAR could either be run : 
- on the matrix produced by cellRanger, with the `chromVAR_per_sample.R` script. You'll need to specify which motif database you'd like to use (HOCOMOCO or TF2DNA)
- on the `fragments.tsv.gz` file, with the `chromVAR_HSC_BAM.R` script. The advantage to this method is that you can use another set of peaks (differentially accessible sites in a specific cluster, sites identified as disrupted by Cicero, ...). I'm using here the `getCountsFromFrags` function from `utils/tools_factory.R` that mimics a custom function written by Caleb Lareau implementation (https://github.com/caleblareau). 

#### Raw Processing 

Contains basic functions to subset the BAM files by a list of barcodes (`subsetBAM.sh`), generate BigWig files for IGV visualisation (`bamnormalizebw_loop.bw`), indexing BAM files and peak calling using MACS2. 

#### HOMER

Peak calling was done on the subsetted BAM files by cluster. The bed produced were then analyzed by HOMER for de novo motif finding and known motif finding

### `Analysis`folder

Basically contains the scripts used for clustering the cells and generating the heatmap of the AUC scores (`big_analysis.R`), to look at the differential accessibility of the motifs using chromVAR (`chromvar_analysis.R`) or HOMER (`deNovoHomer.R` or `knownHomerAnalysis`). 

### `utils` folder

contains the `tools_factory.R`script that contains a few helper functions. I usually load it at the beginning of a lot of my analysis scripts. 

### `ImmGen Preparation`folder







