library(readr)
library(ggplot2)
library(ggrepel)
library(colorblindr)
library(ggpubr)
library(gridExtra)
library(dplyr)

colnames = c('Consensus',
             'p.value',
             'log.p.value',
             'q.value',
             'nb.seq.w.motif',
             'percent.seq.w.motif',
             'bg.nb.seq.w.motif',
             'bg.percent.seq.w.motif')

##### Importing the results for TET2 peaks

tet_results = read.delim('~/Documents/LandauLab/mount/Projects/DNAme/analysis/direct_analysis/HSCs/HOMER/enhancer/enhancer_HCMC_0.0001.TET2.allBG_CORE/knownResults.txt', sep='\t')
colnames(tet_results) = c('MotifName', paste0('TET_', colnames))
tet_results$TET_percent.seq.w.motif = as.numeric(gsub(tet_results$TET_percent.seq.w.motif, pattern='%', replacement=''))
tet_results$TET_bg.percent.seq.w.motif = as.numeric(gsub(tet_results$TET_bg.percent.seq.w.motif, pattern='%', replacement=''))
tet_results$TET_score = tet_results$TET_percent.seq.w.motif/tet_results$TET_bg.percent.seq.w.motif
tet_results$TET_Rank = 1:nrow(tet_results)
rownames(tet_results) = tet_results$MotifName

ggplot(tet_results) + geom_point(aes(x=TET_Rank, y=TET_score), col='navy') + theme_classic() + 
  labs(x='HOMER Rank', y='Custom score') + ggtitle('Custom score VS HOMER Ranks for TET2 results')

dnmt_results = read.delim('~/Documents/LandauLab/mount/Projects/DNAme/analysis/direct_analysis/HSCs/HOMER/enhancer/enhancer_HCMC_0.0001.DNMT3A.allBG_CORE/knownResults.txt', sep='\t')
colnames(dnmt_results) = c('MotifName', paste0('DNMT_', colnames))
dnmt_results$DNMT_percent.seq.w.motif = as.numeric(gsub(dnmt_results$DNMT_percent.seq.w.motif, pattern='%', replacement=''))
dnmt_results$DNMT_bg.percent.seq.w.motif = as.numeric(gsub(dnmt_results$DNMT_bg.percent.seq.w.motif, pattern='%', replacement=''))
dnmt_results$DNMT_score = dnmt_results$DNMT_percent.seq.w.motif/dnmt_results$DNMT_bg.percent.seq.w.motif
dnmt_results$DNMT_Rank = 1:nrow(dnmt_results)
rownames(dnmt_results) = dnmt_results$MotifName

ggplot(dnmt_results) + geom_point(aes(x=DNMT_Rank, y=DNMT_score), col='navy') + theme_classic() + 
  labs(x='HOMER Rank', y='Custom score') + ggtitle('Custom score VS HOMER Ranks for TET2 results')

# Merging the results into a master table 

dnmt_results.sorted = dnmt_results[as.character(rownames(tet_results)), ]

summary(rownames(tet_results)==rownames(dnmt_results.sorted))

homer_ranks = cbind(tet_results,
                    dnmt_results.sorted[, colnames(dnmt_results.sorted)[3:(length(colnames)+3)]])
homer_ranks$MotifName = gsub(rownames(homer_ranks), pattern='_MOUSE.*', replacement = '')
head(homer_ranks)

# Highlight for specific TFs

favorite_TFs = c('GATA1', 'TAL1', 'DDIT3', 'SPI1', 'CEBPA', 'MYC', 'IRF8', 'FLI1', 'ERG', 'KLF3', 
                 'GFI1B', 'RUNX2', 'SOX4', 'MEIS1', 'MEF2C', 'NFIA', 'KLF1', 'GATA2')
homer_ranks$Favorite = homer_ranks$MotifName %in% favorite_TFs

dnmt_score = ggplot() + 
  geom_point(data=homer_ranks, 
             aes(x=DNMT_Rank, y=DNMT_score, fill=Favorite), pch=21) +
  scale_fill_manual(values=c('navy', 'yellow')) +
  geom_text_repel(data=homer_ranks, 
                  aes(x=DNMT_Rank, y=DNMT_score, label=ifelse(Favorite==TRUE, MotifName, ''))) + 
  theme_classic() + 
  geom_point(data=homer_ranks[homer_ranks$DNMT_score<=1, ], 
             aes(x=DNMT_Rank, y=DNMT_score), fill='grey', pch=21) +
  geom_hline(yintercept = 1, linetype='dashed') + 
  labs(x='HOMER Rank', y='Custom score') + 
  ggtitle('Custom score VS HOMER Ranks for DNMT3A results')
tet_score = ggplot() + 
  geom_point(data=homer_ranks,
             aes(x=TET_Rank, y=TET_score, fill=Favorite), pch=21) +
  scale_fill_manual(values=c('navy', 'yellow')) +
  geom_text_repel(data=homer_ranks,aes(x=TET_Rank, y=TET_score, label=ifelse(Favorite==TRUE, MotifName, ''))) + 
  theme_classic() + 
  geom_point(data=homer_ranks[homer_ranks$TET_score<=1, ], 
             aes(x=TET_Rank, y=TET_score), fill='grey', pch=21) +
  geom_hline(yintercept = 1, linetype='dashed') + 
  labs(x='HOMER Rank', y='Custom score') + 
  ggtitle('Custom score VS HOMER Ranks for TET2 results')

grid.arrange(dnmt_score, tet_score, ncol=2)

# Differential scoring

ery_tfs = data.frame(Name = toupper(c('Ddit3', 'Gata1', 'Gata2', 'Gfi1b', 'Klf3','Klf1', 'Tal1', 'Nfia')),
                     Bias = 'Ery')

mono_tfs = data.frame(Name=toupper(c('Irf1', 'Irf8', 'SOX4', 'Spi1', 'Runx2', 'Mef2c', 'CEBPA', 'Meis1', 'FLi1', 'Erg')),
                      Bias = 'Mono')

tfs = rbind(ery_tfs, mono_tfs)
tfs = tfs[order(tfs$Name), ]
rownames(homer_ranks) = homer_ranks$MotifName
homer_ranks_2 = homer_ranks[as.character(tfs$Name),]
summary(tfs$Name==homer_ranks_2$MotifName)

new_cols = c('TET_Rank', 'DNMT_Rank', 'TET_score', 'DNMT_score', 'TET_log.p.value', 'DNMT_log.p.value')
tfs[, new_cols] = homer_ranks_2[ ,new_cols]

tfs_ = tfs[(tfs$TET_score>1)&((tfs$DNMT_score>1)), ]

thresholds = c(100, 125, 150, 200, 250, 300)

plots = list()
for (th in thresholds){
  
  tf_tmp = tfs[(tfs$TET_Rank<th)&((tfs$DNMT_Rank<th)), ]
  tf_tmp$Diff_Ranking = scale(log(tf_tmp$TET_score/tf_tmp$DNMT_score))
  
  plots[[paste0('Thresh_', th)]] = ggplot(tf_tmp) + geom_boxplot(aes(x=Bias, y=Diff_Ranking, fill=Bias), alpha=0.9, width=0.25) + theme_classic() +
    geom_point(aes(x=Bias, y=Diff_Ranking), pch=21, size=1) +
    scale_fill_manual(values=palette_OkabeIto[c(5,6)],
                      breaks=c('Ery', 'Mono'),
                      labels=c(paste0('Ery n=', length(tf_tmp[tf_tmp$Bias=='Ery', 'Bias'])),
                               paste0('Mono n=', length(tf_tmp[tf_tmp$Bias=='Mono', 'Bias'])))) +   
    geom_text_repel(aes(x=Bias, y=Diff_Ranking, label=Name), nudge_x = 0.5) + 
    labs(y='Custom differential score') + geom_hline(yintercept = 0, linetype='dashed', alpha=0.9) + ggtitle(paste0('Rank threshold = ', th)) + 
    stat_compare_means(data=tf_tmp, aes(x=Bias, y=Diff_Ranking), method='t.test', label.y=2.5, label.x=1.4, method.args=list(alternative='greater'))
}

do.call(grid.arrange, c(plots, ncol=3))

tfs_$Diff_Ranking = scale(log(tfs_$TET_score/tfs_$DNMT_score))

ggplot(tfs_) + geom_boxplot(aes(x=Bias, y=Diff_Ranking, fill=Bias), alpha=0.9, width=0.25) + theme_classic() +
  geom_point(aes(x=Bias, y=Diff_Ranking), pch=21, size=1) +
  scale_fill_manual(values=palette_OkabeIto[c(5,6)],
                     breaks=c('Ery', 'Mono'),
                     labels=c(paste0('Ery n=', length(tfs_[tfs_$Bias=='Ery', 'Bias'])),
                              paste0('Mono n=', length(tfs_[tfs_$Bias=='Mono', 'Bias'])))) +   
  geom_text_repel(aes(x=Bias, y=Diff_Ranking, label=Name), nudge_x = 0.5) + 
  labs(y='Custom differential score') + geom_hline(yintercept = 0, linetype='dashed', alpha=0.9) + 
  stat_compare_means(data=tfs_, aes(x=Bias, y=Diff_Ranking), method='t.test', label.y=2.5, label.x=1.4, method.args=list(alternative='less'))

t.test(tf_tmp[tf_tmp$Bias=='Ery', 'Diff_Ranking'], 
       tf_tmp[tf_tmp$Bias=='Mono', 'Diff_Ranking'], 'less')


