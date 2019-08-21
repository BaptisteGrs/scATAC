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

source('Projects/DNAme/Scripts/tools_factory.R')

peaks = resize(getPeaks('Projects/DNAme/Cicero/bed_files/disrupted_sites.bed', sort_peaks=TRUE), 200)
fragment_path = 'Projects/DNAme/data/cellRanger/fragments.tsv.gz'
cellinfo = readRDS('Projects/DNAme/analysis/direct_analysis/cellinfo.rds')

fragment_counts = getCountsFromFrags(fragment_file_XX=fragment_path,
                                     peaks_gr=peaks,
                                     barcodes=rownames(cellinfo[(cellinfo$sampleID %in% c(1,2,3,4,5,6))&(cellinfo$type%in%c('HSC')), ]))

fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Mmusculus.UCSC.mm10)

counts_filtered <- filterPeaks(fragment_counts, non_overlapping = TRUE)

path_to_motif = 'Projects/DNAme/data/HOCOMOCO_mouse_favorite_TFs/'

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

print('Match motifs')
motif_ix = matchMotifs(out, counts_filtered, genome=BSgenome.Mmusculus.UCSC.mm10)

dev = computeDeviations(object = counts_filtered, annotations = motif_ix)

saveRDS(dev, paste0('deviations_HSC_cicero.rds'))



