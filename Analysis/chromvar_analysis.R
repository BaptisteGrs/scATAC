library(chromVAR)
library(ggplot2)
library(viridis)
library(FastKNN)
library(colorblindr)
library(ggpubr)
library(pheatmap)

source('~/Documents/LandauLab/mount/Projects/DNAme/Scripts/tools_factory.R')

# Importing the results 
chromVAR_dev = readRDS('~/Documents/LandauLab/mount/Projects/DNAme/analysis/direct_analysis/chromVAR/strat_replicates_1/deviations.Rda')
chromVAR_z_scores = deviationScores(chromVAR_dev) # Get the z-scores

# Importing the cell info (clustering, genotype, sample ID, ...)
cellinfo = readRDS('~/Documents/LandauLab/mount/Projects/DNAme/analysis/direct_analysis/cellinfo.rds')
cellinfo_ds = cellinfo[rownames(cellinfo) %in% colnames(chromVAR_z_scores), ]
chromVAR_z_scores = chromVAR_z_scores[, colnames(chromVAR_z_scores) %in% rownames(cellinfo_ds)]

cellinfo_ds = cellinfo_ds[colnames(chromVAR_z_scores), ]
summary(rownames(cellinfo_ds) == colnames(chromVAR_z_scores))

cellinfo_ds[, rownames(chromVAR_z_scores)] = t(chromVAR_z_scores)
cellinfo_ds$genotype = factor(cellinfo_ds$genotype, 
                              levels=c('DNMT3A', 'WT', 'TET2'))
tfs = rownames(chromVAR_z_scores)

tf = 'GATA1'
ggplot(cellinfo_ds) + geom_boxplot(aes_string(x='type', y=tf, fill='type')) + theme_classic() + labs(y=paste0(tf, ' z-scores'))

favorite_tfs = c('GATA1', 'GATA2', 'TAL1', 'KLF3', 'DDIT3', 'KLF1', 'NFIA', 
                 'SPI1', 'IRF8', 'CEBPA', 'IRF1', 'SOX4', 'MEIS1', 'FLI1', 'ERG')


pheatmap(cor(cellinfo_ds[, favorite_tfs], method='pearson'), viridis(22))

cellinfo_ds$ery_score_1 = cellinfo_ds$GATA1 + cellinfo_ds$TAL1
cellinfo_ds$ery_score_2 = cellinfo_ds$KLF3 + cellinfo_ds$KLF1

cellinfo_ds$mono_score_1 = cellinfo_ds$IRF8 + cellinfo_ds$IRF1 + cellinfo_ds$MEIS1
cellinfo_ds$mono_score_2 = cellinfo_ds$FLI1 + cellinfo_ds$ERG 

cellinfo_ds$ery_score = cellinfo_ds$ery_score_1 + cellinfo_ds$ery_score_2
cellinfo_ds$mono_score = cellinfo_ds$mono_score_1 + cellinfo_ds$mono_score_2

scores = c('ery_score', 'ery_score_1', 'ery_score_2',
           'mono_score_1', 'mono_score_2', 'mono_score')

box.plots = list()
for (tf in c(favorite_tfs)), scores)){
  box.plots[[tf]] = ggplot(cellinfo_ds[cellinfo_ds$type%in%c('HSC', 'MPP'), ]) + 
    geom_jitter(aes_string(x='genotype', y=tf, fill='genotype'), width=0.3, alpha=0.5, pch=21) + 
    geom_boxplot(aes_string(x='genotype', y=tf, fill='genotype'), alpha=0.9) +
    scale_fill_manual(values=palette_OkabeIto[c(8,5)]) + theme_classic() + 
    stat_compare_means(aes_string(x='genotype', y=tf), ref.group = 'WT',
                       label.y = max(cellinfo_ds[cellinfo_ds$type%in%c('HSC'), tf]))
}

gridExtra::grid.arrange(box.plots$KLF3, 
                        box.plots$KLF1, ncol=2)

gridExtra::grid.arrange(box.plots$GATA1,
                        box.plots$GATA2, 
                        box.plots$TAL1,
                        box.plots$SOX4, ncol=2)

gridExtra::grid.arrange(box.plots$IRF1, 
                        box.plots$IRF8,
                        box.plots$MEIS1,
                        box.plots$NFIA, ncol=2)

gridExtra::grid.arrange(box.plots$SPI1,
                        box.plots$ERG, ncol=2)

gridExtra::grid.arrange(box.plots$FLI1,
                        box.plots$CEBPA, 
                        box.plots$DDIT3, ncol=3)

gridExtra::grid.arrange(box.plots$ery_score_1,
                        box.plots$ery_score_2, 
                        box.plots$mono_score_1,
                        box.plots$mono_score_2, 
                        ncol=2)

gridExtra::grid.arrange(box.plots$ery_score,
                        box.plots$mono_score, ncol=2)
