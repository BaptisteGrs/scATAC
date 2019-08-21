library(cisTopic)
library(readr)
library(AUCell)
library(doRNG)
library(plyr)

source('~/Projects/DNAme/Scripts/tools_factory.R')

args = commandArgs(TRUE)
model = loadModel(args[1], args[2]) #path to cisTopic model, nb of topics

pred.matrix = predictiveDistribution(model)

# Obtain signatures
path_to_signatures = paste0(args[3])
Bulk_ATAC_signatures = paste(path_to_signatures,
                             list.files(path_to_signatures),
                              sep='')
labels  <- gsub('.bed', '', list.files(path_to_signatures))
model <- getSignaturesRegions(model, Bulk_ATAC_signatures, labels=labels, minOverlap = 0.4)

# To only keep unique peaks per signature
#model@signatures <- llply(1:length(model@signatures), function (i) model@signatures[[i]][-which(model@signatures[[i]] %in% unlist(as.vector(model@signatures[-i])))]) 
names(model@signatures) <- labels

# Compute cell rankings 
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

model <- MySignatureCellEnrichment(model, aucellRankings, selected.signatures='all', aucMaxRank = 0.1*nrow(aucellRankings), plot=FALSE)


saveRDS(model, args[4])
