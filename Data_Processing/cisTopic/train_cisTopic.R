library(cisTopic)

# Loading helper functions (used here : createcistopic)
source('Projects/DNAme/data/tools_factory.R')

# Loading filtered barcodes 
load('Projects/DNAme/data/filtered.cells_v2.Rdata')

cistopic = createcistopic(path_to_XX='Projects/DNAme/data/cellRanger/filtered_peak_bc_matrix/',
                          sample_list=as.character(filtered.cells),
                          min_cells=100, 
                          name='Hematopoeisis', 
                          save_path='Projects/DNAme/analysis/direct_analysis/cistopic_blank.rds')

cistopic = runModels(cistopic,
                     topic=c(30,40,50,60,80,100),
                     seed=123,
                     nCores=6,
                     burnin=400,
                     iterations=800)

saveRDS(cistopic, 'Projects/DNAme/analysis/direct_analysis/cistopic_w_models.rds')
