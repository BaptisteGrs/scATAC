print('Importing libraries')

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BiocParallel)

register(MulticoreParam(4))


counts <- Read10X_h5('Projects/DNAme/data/cellRanger/filtered_peak_bc_matrix.h5')
metadata <- read.csv(
  file = 'Projects/DNAme/data/cellRanger/singlecell.csv',
  header = TRUE)

hemato <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)

fragment.path <- 'Projects/DNAme/data/cellRanger/fragments.tsv.gz'

hemato <- SetFragments(
  object = hemato,
  file = fragment.path
)

load('Projects/DNAme/data/filtered.cells_v2.Rdata')

filter_cells = as.character(filtered.cells)

hemato <- subset(hemato, 
                 cells=filter_cells)

hemato <- RunTFIDF(hemato)
hemato <- FindTopFeatures(hemato, min.cutoff = 'q0')
hemato <- RunSVD(
  object = hemato,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

print('Neighbors')
hemato <- FindNeighbors(
  object = hemato,
  reduction = 'lsi',
  dims = 1:30
)
print('Finding clusters')
hemato <- FindClusters(
  object = hemato,
  verbose = FALSE
)

print('Extracting gene coordinates')
#extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object 
gene.coords <- genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

print('Building matrix')
# build a gene by cell matrix
gene.activities <- FeatureMatrix(
  fragments = fragment.path,
  features = genebodyandpromoter.coords,
  cells = colnames(hemato),
  chunk = 10
)

print('Convert rownames')
# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- make.unique(gene.key[rownames(gene.activities)])
gene.activities <- gene.activities[rownames(gene.activities)!="",]

print('Add the gene activity matrix to the Seurat object')
#Add the gene activity matrix to the Seurat object as a new assay, and normalize it
hemato[['RNA']] <- CreateAssayObject(counts = gene.activities)
hemato <- NormalizeData(
  object = hemato,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(hemato$nCount_RNA)
)

print('Saving')
saveRDS(hemato, 'Projects/DNAme/analysis/direct_analysis/Signac/hematoSeurat2.rds')

print('Importing')
hemato = readRDS('Projects/DNAme/analysis/direct_analysis/Signac/hematoSeurat2.rds')
DefaultAssay(hemato) <- 'RNA'
load('Projects/DNAme/analysis/direct_analysis/Signac/seurat_rna.Rdata')

print('Transfering...')
transfer.anchors <- FindTransferAnchors(
  reference = FI.integrated,
  query = hemato,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = FI.integrated$clusterName,
  weight.reduction = hemato[['lsi']]
)

hemato <- AddMetaData(object = hemato, metadata = predicted.labels)

print('Transfer done!!')

saveRDS(hemato, 'Projects/DNAme/analysis/direct_analysis/Signac/hematoSeurat_2_ wRNA.rds')


