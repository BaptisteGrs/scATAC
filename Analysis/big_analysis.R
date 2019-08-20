library(cisTopic)
library(gridExtra)
library(viridis)
library(ggplot2)
library(FastKNN)
library(reshape2)
library(igraph)
library(Rtsne)
library(chromVAR)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggrepel)
library(MASS)
library(rstudioapi)
library(readr)
library(colorblindr)

set.seed(123)

source('Documents/LandauLab/mount/Projects/DNAme/Scripts/tools_factory.R')

cis = loadModel('Documents/LandauLab/mount/Projects/DNAme/analysis/direct_analysis/model_signed.new.rds', 50)
cellinfo = cis@cell.data
cellinfo$sampleID = as.factor(strtoi(sapply(strsplit(as.character(rownames(cellinfo)), ''), tail, 1)))
cellinfo$genotype = factor(ifelse(cellinfo$sampleID %in% 1:2, 'DNMT3A', ifelse(cellinfo$sampleID %in% 3:4, 'TET2', 'WT')), 
                           levels=c('DNMT3A', 'WT', 'TET2'))

ggplot(cellinfo) + 
   geom_boxplot(aes(x=sampleID, y=log(nCounts, 10), fill=genotype)) + theme_classic() + 
   scale_fill_manual(values = c('cornflowerblue', 'grey', 'firebrick2')) + 
   labs(y='log10( # reads in peaks )')

ggplot(cellinfo) + 
   geom_bar(aes(x=sampleID, fill=genotype), col='black') + theme_classic() + 
   scale_fill_manual(values = c('cornflowerblue', 'grey', 'firebrick2')) + 
   labs(y='Nb of unique barcodes')

################################################################################################################################
## Clustering 
################################################################################################################################

cell_topic_distr = modelMatSelection(cis, target = 'cell', method = 'Z-score') 
dist.matrix = 1 - cor(cell_topic_distr, method='spearman')
big.neigh = find_k_neighbors(k=20, dist.matrix=dist.matrix)
big.neigh.m <- melt(big.neigh)
big.neigh.m$Var2 <- NULL
colnames(big.neigh.m) = c('from', 'to')
big.graph = graph_from_data_frame(big.neigh.m, directed=F)
big.communities = cluster_louvain(big.graph)
cellinfo$ClusterID = as.factor(big.communities$membership)

################################################################################################################################
## Dimensionality reduction with UMAP
################################################################################################################################

library(umap)
custom.config = umap.defaults
custom.config$random_state = 123

cis.umap = umap(t(cell_topic_distr))
cellinfo[, c('UMAP1', 'UMAP2')] = cis.umap$layout

ggplot(cellinfo) + 
   geom_bar(aes(x=ClusterID, fill=genotype, alpha=as.numeric(sampleID)%%2==0), col='black') + theme_classic() + 
   scale_fill_manual(values = c('cornflowerblue', 'grey', 'firebrick2')) + 
   scale_alpha_discrete(name='Replicate', breaks = c(FALSE, TRUE), labels = c(1, 2), 
                        range=c(0.7, 1)) + 
   labs(y='Nb of cells', x='Cluster ID')

ggplot(cellinfo) + 
   geom_bar(aes(x=ClusterID, fill=genotype, alpha=as.numeric(sampleID)%%2==0), position='fill', col='black') + theme_classic() + 
   scale_fill_manual(values = c('cornflowerblue', 'grey', 'firebrick2')) + 
   scale_alpha_discrete(name='Replicate', breaks = c(FALSE, TRUE), labels = c(1, 2), 
                        range=c(0.7, 1)) + 
   labs(y='Nb of cells', x='Cluster ID')

ggplot(cellinfo) + 
   geom_boxplot(aes(x=ClusterID, y=nCounts, fill=ClusterID), notch=TRUE) + theme_classic() + 
   labs(y='Nb of reads in peaks', x='Cluster ID')

mean_cluster = cellinfo %>% 
                  group_by(ClusterID) %>%
                  summarize(mean_UMAP1 = median(UMAP1), mean_UMAP2 = median(UMAP2)) %>%
                  data.frame() 

ggplot(cellinfo) + 
   geom_point(aes(x=UMAP1, y=UMAP2, col=ClusterID), size=1) + 
   theme_classic() + geom_label(data=mean_cluster, aes(x=mean_UMAP1, y=mean_UMAP2, label=ClusterID), size=5)

################################################################################################################################
# ImmGen signatures
################################################################################################################################

cell_types = read.csv('Documents/LandauLab/mount/Projects/DNAme/data/ImmGen/cell_types.csv', sep=';')
cellinfo[, paste0('Scaled_', c(colnames(cell_types), 'Erythroid', 'Megakaryocyte'))] = scale(cellinfo[, c(colnames(cell_types), 'Erythroid', 'Megakaryocyte')])

cluster_annotations = c('CMP', 'LMPP', 'GMP', 'HSC', 'MPP', 'MEP', 'CLP')
names(cluster_annotations) = 1:7

cellinfo$type = factor(cluster_annotations[as.character(cellinfo$ClusterID)], levels=c('HSC', 'MPP', 'LMPP', 'CLP', 'CMP', 'GMP', 'MEP'))

cellinfo.sorted = cellinfo[order(cellinfo$type), ]
imm.scores.mat.sorted = t(as.matrix(cellinfo.sorted[, paste0('Scaled_', c(colnames(cell_types), 'Erythroid', 'Megakaryocyte'))]))

rownames(imm.scores.mat.sorted) = c(colnames(cell_types), 'Erythroid', 'Megakaryocyte')

cluster_assignment = data.frame("Clusters" = as.character(cellinfo.sorted$type),
                                row.names = rownames(cellinfo.sorted))

hm_colors = colorRampPalette(brewer.pal(7,"Dark2"))(7)
names(hm_colors) = cluster_annotations

cluster_assignment$Clusters = factor(cluster_assignment$Clusters, levels=c('HSC', 'MPP', 'LMPP', 'CLP', 'CMP', 'GMP', 'MEP'))

cluster_assignment = cluster_assignment[order(cluster_assignment$Clusters), , drop=FALSE]

annotation = HeatmapAnnotation(df=cluster_assignment, show_legend = FALSE, border=TRUE,
                               col = list(Clusters=hm_colors),
                               which = 'column', annotation_legend_param = gpar(fontsize=3))

imm.scores.mat.sorted = imm.scores.mat.sorted[, rownames(cluster_assignment[order(cluster_assignment$Clusters), , drop=FALSE])]

#pdf('./Heatmap_AUC.pdf')
Heatmap(scale(imm.scores.mat.sorted),
        cluster_rows = FALSE,
        col=viridis_pal()(22),
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        top_annotation = annotation, column_split = cluster_assignment, column_gap = unit(4, "mm"))
#dev.off()

saveRDS(cellinfo, 'Documents/LandauLab/mount/Projects/DNAme/direct_analysis/cellinfo.rds')