library(chromVAR)
library(chromVARmotifs)
library(JASPAR2018)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
library(SummarizedExperiment)
library(data.table)
library(readr)

register(MulticoreParam(4))
source('~/Projects/DNAme/Scripts/tools_factory.R')

args = (commandArgs(TRUE))
sample = eval(parse(text=args[1]))
name = args[2]
method=args[3]

cellinfo = readRDS('Projects/DNAme/TET2/cellinfo.rds')
reads_matrix = readRDS('Projects/DNAme/TET2/reads_matrix.rds')

# peaks to keep
peak_data = transpose(as.data.table(strsplit(rownames(reads_matrix), '[:-]')))
colnames(peak_data) <- c('chr', 'start', 'end')
atac_range = makeGRangesFromDataFrame(peak_data)
print(length(atac_range))

# cells to keep 
cellinfo = cellinfo[cellinfo$sampleID %in% sample, ]
cell_to_keep = colnames(reads_matrix) %in% rownames(cellinfo)

# Filtering the matrix
matrix = reads_matrix[, cell_to_keep]
print(dim(matrix))


colData = cellinfo[colnames(matrix), ]
colData$depth = Matrix::colSums(matrix)

print('Initializing chromVAR objects...')
chrom_counts <- SummarizedExperiment(assays=list(counts=matrix),
                                    rowRanges=atac_range, colData=colData)

print('Filtering out peaks with 0 cells covered')
chrom_counts <- filterPeaks(chrom_counts) 

chrom_counts <- addGCBias(chrom_counts, genome=BSgenome.Mmusculus.UCSC.mm10)

if (method=='HOCOMOCO'){
   path_to_motif = 'Projects/DNAme/data/HOCOMOCO_mouse/'

   file.names <- dir(path_to_motif, pattern =".motif|.mat|.txt")
   myMotifs = list()
   for(i in 1:length(file.names)){
         myMotifs[[i]] <- as.data.frame(read_delim(paste0(path_to_motif,file.names[i]), 
                                                   "\t"))
         colnames(myMotifs[[i]]) =c("A", "C", "G", "T")
         rownames(myMotifs[[i]]) = c(1:nrow(myMotifs[[i]]))
   }
   names(myMotifs) = gsub(file.names, pattern = ".motif|.mat|.txt", replacement = "")
   names(myMotifs) = toupper(names(myMotifs))

   to_PWMatrix = function(liste){
      out = PWMatrixList()
      for (i in 1:length(liste)){
         mat = t(as.matrix(liste[[i]]))*100000
         name = names(myMotifs)[i]
         mode(mat) = 'integer'
         out[[name]] = PWMatrix(ID = name, name = names(myMotifs)[i], profileMatrix = toPWM(mat, type = 'prob'))
      }
      return(out)
   }

   out = to_PWMatrix(myMotifs)
} else if (method=='TF2DNA'){
   path_to_motif = 'Projects/DNAme/data/TF2DNA_subset/'
   file.names <- dir(path_to_motif, pattern =".motif|.mat|.txt")
   myMotifs = list()
   for(i in 1:length(file.names)){
      myMotifs[[i]] <- as.data.frame(read_delim(paste0(path_to_motif,file.names[i]), 
                                                 " ", escape_double = FALSE, col_names = FALSE, 
                                                 trim_ws = TRUE))
      myMotifs[[i]]$X1 = NULL
      myMotifs[[i]] = t(as.matrix(myMotifs[[i]]))
      colnames(myMotifs[[i]]) =c("A", "C", "G", "T")
      rownames(myMotifs[[i]]) = c(1:nrow(myMotifs[[i]]))
   }
   names(myMotifs) = gsub(file.names, pattern = ".motif|.mat|.txt", replacement = "")
   names(myMotifs) = toupper(names(myMotifs))

   to_PWMatrix = function(liste){
      out = PWMatrixList()
      for (i in 1:length(liste)){
         mat = t(as.matrix(liste[[i]]))*100000
         name = names(myMotifs)[i]
         mode(mat) = 'integer'
         out[[name]] = PWMatrix(ID = name, name = names(myMotifs)[i], profileMatrix = toPWM(mat, type = 'prob'))
      }
      return(out)
   }

   out = to_PWMatrix(myMotifs)
}


print('Match motifs')
motif_ix = matchMotifs(out, chrom_counts, genome=BSgenome.Mmusculus.UCSC.mm10)

dev = computeDeviations(object = chrom_counts, annotations = motif_ix)
saveRDS(dev, paste0('deviations_', name,'.Rda'))




