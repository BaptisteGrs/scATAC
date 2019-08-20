library(rstudioapi)
library(igraph)
library(tibble)
library(GenomicRanges)

# Setting the working directory 
WD = dirname(getSourceEditorContext()$path)
setwd(WD)

##--------------
# ImmGen regions
##--------------

immgen.table = read.csv('./ImmGenATAC18_AllOCRsInfo.csv') # raw table from ImmGen
head(immgen.table)

# Expanding 250 bp upstream and downstream
immgen.table = add_column(immgen.table, start=immgen.table$Summit-250, .after = 'Summit')
immgen.table = add_column(immgen.table, end=immgen.table$Summit+250, .after = 'start')

# Finding overlapping peaks
imm_ranges = makeGRangesFromDataFrame(immgen.table, keep.extra.columns = TRUE)
overlaps = findOverlaps(imm_ranges)

ind1 = queryHits(overlaps)
ind2 = subjectHits(overlaps)
over_indices = as.data.frame(table(ind1))

immgen.table = add_column(immgen.table, Freq = over_indices$Freq, .after='end')
immgen.table = add_column(immgen.table, ID=seq.int(nrow(immgen.table)), .after='ImmGenATAC1219.peakID')

graph = readRDS('./peaks_graph.rds') ##### job run on cluster to eliminate the replicates
vertices = strtoi(V(graph)$name)

immgen.table.final = immgen.table[(immgen.table$Freq==1)|(immgen.table$ID %in% vertices), ]
colnames(immgen.table.final)[3] = 'chr'

imm_ranges_final = makeGRangesFromDataFrame(immgen.table.final)

length(queryHits(findOverlaps(imm_ranges_final)))==dim(immgen.table.final)[1]
# No more overlaps. 407475  peaks, 500 width

# Fixing the coordinates 

chr_size = read.table('./mm10_sizes.txt')
colnames(chr_size) = c('chr', 'max_size')
immgen.table.final[which(immgen.table.final$start <= 0), 'start'] = 1

for (row in 1:dim(immgen.table.final)[1]){
  if (immgen.table.final[row, 'end'] > chr_size$max[chr_size[, 'chr']==as.character(immgen.table.final[row, 'chr'])]){
    immgen.table.final[row, 'end'] = chr_size$max[chr_size[, 'chr']==as.character(immgen.table.final[row, 'chr'])]
  }
}

##------------------
# Merging cell types
##------------------

cell_types = read.csv('./cell_types.csv', sep=';')

summarized_immgen = data.frame(matrix(ncol=dim(cell_types)[2]+3, nrow=dim(immgen.table.final)[1]))
colnames(summarized_immgen) = c('chr', 'start', 'end', colnames(cell_types))
summarized_immgen[, c('chr', 'start', 'end')] = immgen.table.final[, c('chr', 'start', 'end')]
for (celltype in colnames(cell_types)){
  subtypes = levels(cell_types[, celltype])
  tmp = immgen.table.final[, subtypes, drop=FALSE]
  summarized_immgen[, celltype] = rowMeans(tmp)
}

pathToSave = './BulkProfiles/Top30k/'
top_n = 30000 # top 60000 regions

# filter out chromosomes that are not in the cellRanger output 
peaks_cellranger = read.table('../../analysis/filtered_peaks.bed', sep='\t')
summarized_immgen = summarized_immgen[summarized_immgen$chr %in% levels(unique(peaks_cellranger$V1)), ]
dim(summarized_immgen)

saveRDS(summarized_immgen, './summarized_immgen.rds')
summarized_immgen = readRDS('./summarized_immgen.rds')

for (celltype in colnames(cell_types)){
  print(celltype)
  tmp = summarized_immgen[order(-summarized_immgen[, celltype]), ]
  write.table(x = tmp[1:top_n, c('chr', 'start', 'end', celltype)], file=paste0(pathToSave, celltype, '.bed'),
              quote=FALSE, row.names = FALSE,
              col.names=FALSE)
}

##--------------
# ENCODE regions
##--------------

# For erythroid and megakaryocyte profiles 
library(rtracklayer)

extraCols_narrowPeak = c(signalValue = "numeric",
                          pValue = "numeric",
                          qValue = "numeric",
                          peak = "integer")

erythroid_range = sort(import('./ENCODE/ENCFF005XHN.bed',format = "BED", extraCols = extraCols_narrowPeak))
megakaryocyte_range = sort(import('./ENCODE/ENCFF066SZX.bed', format = "BED", extraCols = extraCols_narrowPeak))

### Filter the one with best p-value
ery.df = as.data.frame(erythroid_range)
ery.df = ery.df[order(-ery.df$pValue), ]
ery_range = makeGRangesFromDataFrame(ery.df, keep.extra.columns = TRUE)
ery_range = unique(ery_range)
ery.df.export = as.data.frame(ery_range)

mega.df = as.data.frame(megakaryocyte_range)
mega.df = mega.df[order(-mega.df$pValue), ]
mega_range = makeGRangesFromDataFrame(mega.df, keep.extra.columns = TRUE)
mega_range = unique(mega_range)
mega.df.export = as.data.frame(mega_range)

### Re-center the peaks and expand 250 bp upstream and downstream
ery.df.export$width = ery.df.export$end - ery.df.export$start
ery.df.export$check = ery.df.export$width/2
ery.df.export$peak_coord = ery.df.export$start + ery.df.export$peak
ery.df.export$start = ery.df.export$peak_coord -250
ery.df.export$end = ery.df.export$peak_coord + 250
for (row in 1:dim(ery.df.export)[1]){
  if (ery.df.export[row, 'end'] > chr_size$max[chr_size[, 'chr']==as.character(ery.df.export[row, 'seqnames'])]){
    ery.df.export[row, 'end'] = chr_size$max[chr_size[, 'chr']==as.character(ery.df.export[row, 'seqnames'])]
  }
}

mega.df.export$width = mega.df.export$end - mega.df.export$start
mega.df.export$check = mega.df.export$width/2
mega.df.export$peak_coord = mega.df.export$start + mega.df.export$peak
mega.df.export$start = mega.df.export$peak_coord -250
mega.df.export$end = mega.df.export$peak_coord + 250
for (row in 1:dim(mega.df.export)[1]){
  if (mega.df.export[row, 'end'] > chr_size$max[chr_size[, 'chr']==as.character(mega.df.export[row, 'seqnames'])]){
    mega.df.export[row, 'end'] = chr_size$max[chr_size[, 'chr']==as.character(mega.df.export[row, 'seqnames'])]
  }
}

summary(ery.df.export)
summary(mega.df.export)

ery.df.export = ery.df.export[order(-ery.df.export$score), ]
write.table(x = ery.df.export[1:top_n, c('seqnames', 'start', 'end', 'score')], file=paste0(pathToSave, 'Erythroid.bed'),
            quote=FALSE, row.names = FALSE,
            col.names=FALSE)

mega.df.export = mega.df.export[order(-mega.df.export$score), ]
write.table(x = mega.df.export[1:top_n, c('seqnames', 'start', 'end', 'score')], file=paste0(pathToSave, 'Megakaryocyte.bed'),
            quote=FALSE, row.names = FALSE,
            col.names=FALSE)


