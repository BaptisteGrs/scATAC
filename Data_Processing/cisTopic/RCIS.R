library(cisTopic)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)

source('Projects/DNAme/Scripts/tools_factory.R')

cis = readRDS('Projects/DNAme/analysis/model_signed_2.Rds')

txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

cis <- annotateRegions(cis, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb='org.Mm.eg.db')
cis <- getRegionsScores(cis, method='NormTop', scale=TRUE)
cis <- binarizecisTopics(cis, thrP=0.975, plot=FALSE)

path_to_feather = 'Projects/DNAme/data/mm9-tss-centered-10kb-10species.mc9nr.feather'

mm10Tomm9.chain  <- import.chain("mm10Tomm9.over.chain")

# Obtain liftOver dictionary (as list)
print('cistopic pipeline')
mm10_coord <- cis@region.ranges
mm10_to_mm9_list <- liftOver(mm10_coord, mm10Tomm9.chain)

cis <- binarizedcisTopicsToCtx(cis, liftOver=mm10_to_mm9_list, genome='mm9')

cis <- scoredRegionsToCtx(cis, liftOver=mm10_to_mm9_list, genome='mm9')

library(RcisTarget)
library(plyr)
topicsList <- cis@binarized.regions.to.Rct
motifRankings <- importRankings(path_to_feather)
ctxregions <- colnames(getRanking(motifRankings))[-1]
topicsList <- llply(1:length(topicsList), function(i) topicsList[[i]][which(topicsList[[i]] %in% ctxregions)])
names(topicsList) <- names(cis@binarized.regions.to.Rct)
saveRDS(topicsList, 'Projects/DNAme/analysis/topicsList.rds')
saveRDS(motifRankings, 'Projects/DNAme/analysis/motifRankings.rds')

cis <- topicsRcisTarget(cis, genome='mm9',
                        path_to_feather, 
                        reduced_database=TRUE,
                        nesThreshold=3,
                        rocthr=0.005,
                        maxRank=20000)
                        #nCores=1)

cis<- getCistromes(cis, annotation = 'Both', nCores=5)

saveRDS(cis, 'Projects/DNAme/analysis/model_rcis.rds')

