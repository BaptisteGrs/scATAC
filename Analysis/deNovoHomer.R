# Pseudo-bulk profiles were generated for the HSC cluster by genotype, producing 3 bam files

# Peak calling was run separately on these 3 bam files, with the following parameters : 
# --broad --broad-cutoff 0.1 -m 2 100 

# Homer motif finding function was used on the 2 KO genotypes bed files produced : -size 200 -mask

library(readr)
library(ggplot2)
library(ggrepel)
library(colorblindr)
library(ggpubr)
library(gridExtra)
library(dplyr)

##########################
#### DE NOVO MOTIF FINDING
##########################

de_novo_tet_path = '~/Documents/LandauLab/mount/Projects/DNAme/analysis/direct_analysis/HSCs/HOMER/TET2_ALLCELLSBG/homerResults/'
tet.motifs = list.files(de_novo_tet_path, pattern='motif...motif|motif..motif')

de_novo_dnmt_path = '~/Documents/LandauLab/mount/Projects/DNAme/analysis/direct_analysis/HSCs/HOMER/DNMT3A_ALLCELLSBG/homerResults/'
dnmt.motifs = list.files(de_novo_dnmt_path, pattern='motif...motif|motif..motif')

read_motif = function(path_to_motif){
  con = file(path_to_motif, open='r')
  firstLine = readLines(con,n=1, warn=FALSE)
  p.value = as.numeric(unlist(strsplit(firstLine, '\t'))[4])
  name = paste(unlist(strsplit(unlist(strsplit(unlist(strsplit(firstLine, '\t'))[2], 'BestGuess:'))[2], '/'))[c(1,2)], collapse='_')
  cpgScore = NULL
  matrix = matrix(data=as.numeric(unlist(strsplit(unlist(strsplit(readLines(con, n = -1, warn = FALSE), '\n')), '\t'))), ncol=4, byrow=TRUE)
  colnames(matrix) = c('A', 'C', 'G', 'T')
  
  for (m in 1:(nrow(matrix)-1)){
    cpgScore = c(cpgScore, matrix[m, 2]*matrix[m+1, 3])
    #cpgScore = c(cpgScore,(1 - (sum(matrix[m,c(1,3,4)]) + matrix[m,2]*sum(matrix[m+1, c(1,2,4)]))))
  }
  close(con)
  return(list(name=name, p.value= p.value, matrix=matrix, cpgScore=cpgScore))
}

df.motifs.tet = data.frame(matrix(nrow=length(tet.motifs), ncol=5))
colnames(df.motifs.tet) = c('Name', 'MaxCpG', 'MeanCpG', 'p.value', 'Rank')

for (i in 1:length(tet.motifs)){
  motif = read_motif(paste0(de_novo_tet_path, tet.motifs[i]))
  df.motifs.tet[i, ]$Name = motif$name
  df.motifs.tet[i, ]$MaxCpG = max(motif$cpgScore)
  df.motifs.tet[i, ]$MeanCpG = mean(motif$cpgScore)
  df.motifs.tet[i, ]$p.value = motif$p.value
  df.motifs.tet[i, ]$Rank = strtoi(gsub(tet.motifs[i], pattern='motif|\\.', replacement = ''))
}

df.motifs.dnmt = data.frame(matrix(nrow=length(dnmt.motifs), ncol=5))
colnames(df.motifs.dnmt) = c('Name', 'MaxCpG', 'MeanCpG', 'p.value', 'Rank')

for (i in 1:length(dnmt.motifs)){
  motif = read_motif(paste0(de_novo_dnmt_path, dnmt.motifs[i]))
  df.motifs.dnmt[i, ]$Name = motif$name
  df.motifs.dnmt[i, ]$MaxCpG = max(motif$cpgScore)
  df.motifs.dnmt[i, ]$MeanCpG = mean(motif$cpgScore)
  df.motifs.dnmt[i, ]$p.value = motif$p.value
  df.motifs.dnmt[i, ]$Rank = strtoi(gsub(dnmt.motifs[i], pattern='motif|\\.', replacement = ''))
}

# Discretize CpG content
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

df.motifs.dnmt$Geno = as.factor('DNMT3A')
df.motifs.tet$Geno = as.factor('TET2')

df.motifs = rbind(df.motifs.dnmt, df.motifs.tet)
df.motifs$Rank = as.numeric(df.motifs$Rank)

df.motifs = df.motifs[df.motifs$p.value<exp()]

summary(df.motifs)

df.motifs = mutate(df.motifs, quantile_rank=ntile(df.motifs$MeanCpG, 4))
df.motifs$quantile_rank = as.factor(df.motifs$quantile_rank)
df.motifs$Name = gsub(df.motifs$Name, pattern='_MOUSE.*', replacement = '')
  
# BOX PLOT 

df.motifs$CpG_Category = factor(ifelse(df.motifs$quantile_rank %in% c(1), 'Low CpG',
                                ifelse(df.motifs$quantile_rank %in% c(4), 'High CpG', 
                                       'Average CpG')), 
                                levels=c('Low CpG', 'Average CpG', 'High CpG'))

ggplot(df.motifs[df.motifs$CpG_Category!='Average CpG', ]) + geom_boxplot(aes(x=CpG_Category, 
                                     y=Rank, 
                                     fill=Geno), 
                                     alpha=0.9)  + 
  geom_jitter(aes(x=CpG_Category, 
                 y=Rank, fill=Geno),position=position_dodge(width=0.75), pch=21) +
  scale_fill_manual(values=palette_OkabeIto[5:6]) + 
  theme_classic() + labs(y='HOMER Rank') +
  ylim(rev(c(-1, max(range(df.motifs$Rank))))) + 
  scale_x_discrete(name='CpG Content', breaks=c('Low CpG', 'High CpG'), labels=c('Low CpG \n (Bottom Quartile)', 'High CpG \n (Top Quartile)'))


fisher.test(x=df.motifs[df.motifs$CpG_Category != 'Average CpG', 'Geno'], 
            y=df.motifs[df.motifs$CpG_Category != 'Average CpG', 'CpG_Category'])

# STAIR PLOT

ggplot(df.motifs) +
  geom_bar(aes(x=quantile_rank, fill=Geno), col='black', position='fill') + 
  scale_fill_manual(name='Genotype', values=palette_OkabeIto[c(5,6)]) +
  labs(x='Mean p(CpG) by quartile', y='# of de novo motifs \n TET2 VS DNMT3A') + 
  scale_x_discrete(breaks=1:10, labels=paste0('Q', 1:10)) +
  theme_classic()

ggplot(df.motifs, aes(MeanCpG, col=Geno)) + 
  stat_ecdf(geom='step') + theme_classic() + labs(y='CDF')

wilcox.test(x=df.motifs[df.motifs$Geno=='DNMT3A', 'MeanCpG'], 
        y=df.motifs[df.motifs$Geno=='TET2', 'MeanCpG'], alternative = 'greater', paired=FALSE)

