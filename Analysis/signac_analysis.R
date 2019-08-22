library(ggplot2)
library(colorblindr)
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)
library(ggplot2)

hemato = readRDS('Documents/LandauLab/mount/Projects/DNAme/analysis/direct_analysis/cellinfo.rds')

seurat_results = readRDS('Documents/LandauLab/mount/Projects/DNAme/analysis/direct_analysis/Signac/metadata.rds')
seurat_umap = readRDS('Documents/LandauLab/mount/Projects/DNAme/analysis/direct_analysis/Signac/umap_signac.rds')

summary(rownames(hemato)==rownames(seurat_umap))

hemato[, colnames(seurat_results)] = seurat_results
hemato[, paste0('Seurat_', colnames(seurat_umap))] = seurat_umap

head(hemato)

cis.umap = ggplot(hemato) + geom_point(aes(x=UMAP1, y=UMAP2, fill=type), pch=21, size=2) + 
  theme_classic() + ggtitle('Dim reduction with cisTopic')

lsi.umap = ggplot(hemato) + geom_point(aes(x=Seurat_UMAP_1, y=Seurat_UMAP_2, fill=type), pch=21, size=2) + 
  theme_classic() + ggtitle('Dim reduction with Seurat (TF-IDF/LSI)')

gridExtra::grid.arrange(cis.umap, lsi.umap, ncol=2)

# Heatmap

hemato$type = factor(hemato$type, levels=c('HSC', 'MPP', 'LMPP', 'CLP', 'CMP', 'GMP', 'MEP'))
hemato.sorted = hemato[order(hemato$type), ]

cols_to_keep = colnames(hemato.sorted)[grep(colnames(hemato.sorted), pattern='prediction.score')]
seurat_score = hemato.sorted[, cols_to_keep, drop=FALSE]

seurat_score$prediction.score.Osteoclast = NULL
seurat_score$prediction.score.Platelet = NULL
seurat_score$prediction.score.max = NULL
colnames(seurat_score) = gsub(colnames(seurat_score), pattern='prediction.score.', replacement = '')
seurat_score = seurat_score[, c("HSC.1"
                              ,"HSC.2"
                              ,"HSC.3"
                              ,"HSC_cc"
                              ,"IMP1"
                              ,"IMP2"
                              ,"Mono.1"
                              ,"Mono.2"
                              ,"Mono.3"
                              ,"Ba.1"
                              ,"Ba.2"
                              ,"Ba.3"
                              ,"Ba.4"
                              ,"Eo"
                              ,"CLP"
                              ,"B.cell"
                              ,"T.cell"
                              ,"T.cell.Cd3d."
                              ,"Ery.1"
                              ,"Ery.2"
                              ,"Ery.3"
                              ,"Ery.4"
                              ,"Neu.1"
                              ,"Neu.2"
                              ,"Neu.3"
                              ,"Neu.4"
                              ,"MkP.1"
                              ,"MkP.2"
                              ,"E.B"
                              ,"NA")
                            ]

seurat_score.m = scale(t(scale(as.matrix(seurat_score))))

cluster_assignment = data.frame("Clusters" = as.character(hemato.sorted$type),
                                row.names = rownames(hemato.sorted))

hm_colors = palette_OkabeIto[1:7]
annotation = HeatmapAnnotation(df=cluster_assignment, show_legend = FALSE, border=TRUE,
                               col = list(Clusters=c("HSC"=hm_colors[1],
                                                     "MPP"=hm_colors[2],
                                                     "LMPP"=hm_colors[3],
                                                     "CLP"=hm_colors[4],
                                                     "CMP"=hm_colors[5],
                                                     "GMP"=hm_colors[6],
                                                     "MEP"=hm_colors[7])),
                               which = 'column', annotation_legend_param = gpar(fontsize=3))

cluster_assignment$Clusters = hemato.sorted$type

ComplexHeatmap::Heatmap(seurat_score.m,
        cluster_rows = FALSE,
        col=viridis_pal()(22),
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        top_annotation = annotation, column_split = cluster_assignment, column_gap = unit(4, "mm"))
        #heatmap_legend_param = list(
         # title = "Scaled Prediction Score", at = -4:4, legend_height = unit(6, "cm")))

df.per = data.frame(matrix(nrow=7, ncol =length(unique(hemato$predicted.id))))
rownames(df.per) = unique(hemato.sorted$type)
colnames(df.per) = colnames(seurat_score)

for (mycluster in rownames(df.per)){
  for (rna_cluster in colnames(df.per)){
    nb = dim(hemato.sorted[(hemato.sorted$type==mycluster)&(gsub((hemato.sorted$predicted.id), pattern='[-|_|/\\+]', replacement='.')==rna_cluster), ])[1]
    df.per[mycluster, rna_cluster] = 100*nb/dim(hemato.sorted[(hemato.sorted$type==mycluster), ])[1]  
  }
}

library(pheatmap)

pheatmap(df.per, cluster_rows = FALSE, cluster_cols = FALSE, color = viridis(20), angle_col = 45)

hsc = hemato[(hemato$predicted.id%in%c('HSC-1', 'HSC-2', 'HSC-3'))&(hemato$type%in%c('HSC')), ]

write.table(rownames(hsc),
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE,
            file='Documents/LandauLab/mount/Projects/DNAme/analysis/direct_analysis/HSCs/chromVAR/barcodes.txt')


ggplot(hemato) + 
  geom_bar(aes(x=sample, fill=sample), col='black') +
  theme_classic() + 
  scale_fill_manual(values=palette_OkabeIto[sort(rep(c(5,6,8), 2))]) + 
  labs(x='Sample ID', y='Nb of cells')

ggplot(hemato) + 
  geom_boxplot(aes(x=sample, y=nAcc, fill=sample), notch=TRUE, col='black') +
  theme_classic() + 
  scale_fill_manual(values=palette_OkabeIto[sort(rep(c(5,6,8), 2))]) + 
  labs(x='Sample ID', y='Nb of accessible sites')

#######@ TET2 Exon 3

cells_tet1 = read.delim('~/Documents/LandauLab/mount/Projects/DNAme/data/cellRanger/DNMT3A_atac_bam/BAM_QC_TET2_EXON3/cells_tet2/reads_barcode_tet2_1',
                        header=FALSE, sep=' ')

cells_tet2 = read.delim('~/Documents/LandauLab/mount/Projects/DNAme/data/cellRanger/DNMT3A_atac_bam/BAM_QC_TET2_EXON3/cells_tet2/reads_barcode_tet2_2',
                        header=FALSE, sep=' ')

colnames(cells_tet1) = c('nb.reads', 'barcode')
colnames(cells_tet2) = c('nb.reads', 'barcode')

cells_tet1$barcode = paste0(cells_tet1$barcode, '-3')
cells_tet2$barcode = paste0(cells_tet2$barcode, '-4')

head(hemato)
hemato$outliers = as.factor(ifelse((rownames(hemato) %in% cells_tet1$barcode)|(rownames(hemato) %in% cells_tet2$barcode), 'Outliers', 'Normal'))

ggplot(hemato) + geom_point(aes(x=UMAP1, y=UMAP2),
                                col=palette_OkabeIto[8], alpha=0.8) + 
  geom_point(data=hemato[hemato$outliers=='Outliers', ], 
             aes(x=UMAP1, y=UMAP2),
             col=palette_OkabeIto[5], size=2) + 
  theme_classic() 


