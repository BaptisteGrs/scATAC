suppressMessages(library(readr))
suppressMessages(library(rtracklayer))
library(tibble)
library(igraph)

table0 <- suppressMessages(read_csv('Projects/DNAme/data/ImmGen/ImmGenATAC18_AllOCRsInfo.csv'))

start <- table0$Summit - 250
end <- table0$Summit + 250

table1 <- add_column(table0, start=start, .after='Summit')
table2 <- add_column(table1, end=end, .after='start')

imm_range0 <- makeGRangesFromDataFrame(table2, keep.extra.columns = T)

overlaps <- findOverlaps(imm_range0, imm_range0)

ind1 <- queryHits(overlaps)
ind2 <- subjectHits(overlaps)
indices <- data.frame(table(ind1))
uni_indices <- subset(indices, indices$Freq > 1)
duplicate_indices <- uni_indices$ind1

dupl_table <- table2[duplicate_indices, ]
dupl_table <- add_column(dupl_table, Freq=uni_indices$Freq, .after='end')

indices_ <- data.frame(ind1, ind2)
indices_ <- subset(indices_, (indices_$ind1!=indices_$ind2))

p_ref_table <- data.frame(indices_$ind1, table2[indices_$ind1, '_-log10_bestPvalue' ])
colnames(p_ref_table) <- c('indx', 'pvalue')
p_ref_table <- unique(p_ref_table[order(-p_ref_table$pvalue), ]) #sorting the peak per p-value

g <- graph_from_data_frame(indices_, directed=F)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

indices_2 <- indices_
p_ref_table_2 <- p_ref_table

for (i in 1:dim(p_ref_table_2)[1]){
    idx <- p_ref_table_2$indx[i]
    if (i%%100==0){
        print(100*i/dim(p_ref_table_2)[1])
    }
    g <- tryCatch(
        {delete.vertices(g, neighbors(graph = g, v = toString(idx)))
         },
        error=function(cond){
        return(g)})
}
g

saveRDS(g, 'Projects/DNAme/data/ImmGen/peaks_graph.Rds')
