### Set of helper functions 

loadModel <- function(model_path, nb_topic=NULL){ 
    # Load cisTopic model from a .Rds file with a certain number of topic

    print('Reading model...')
    model <- readRDS(model_path)
    print('Loaded!')
    if(is.null(nb_topic)){
        model <- selectModel(model)
    }
    else{
        model <- selectModel(model, select=nb_topic)
    }
    colnames(model@cell.data) <- make.names(colnames(model@cell.data), unique=TRUE)
    head(model@cell.data)
    return(model)
}

MySignatureCellEnrichment <- function(object, aucellRankings, selected.signatures='all', nCores = 1, aucMaxRank = 0.03*nrow(aucellRankings), plot=TRUE, ...){
    # There was a typo at the time in the cisTopic source code that was preventing the usage of signatureCellEnrichment function
    # This is just the typo-corrected version

    # Check info
    if((length(object@signatures) < 1)&&(selected.signatures!='annotation')){
        stop('Please, run getSignaturesRegions() first.')
    }

    if (selected.signatures!='annotation'){
        signatures <- object@signatures
        if (is.null(signatures)){
            stop('Please run getSignaturesRegions() first.')
        }
        if(selected.signatures != 'all'){
            if (sum(selected.signatures %in% names(signatures)) != length(selected.signatures)){
                stop('Check whether the selected signatures have been stored in object@signatures.')
            }
            signatures <- signatures[selected.signatures]
        }
    }
    else{
    signatures <- split(object@region.data, object@region.data$annotation)    
    if (is.null(signatures)){
    stop('Please run annotateRegions() first.')
    }
    signatures <- lapply(signatures, function(x) rownames(x))
    }

    modulesAUC <- AUCell_calcAUC(signatures, aucellRankings, nCores=nCores, aucMaxRank=aucMaxRank)
    enrichMatrix <- t(getAUC(modulesAUC))
    rownames(enrichMatrix) <- object@cell.names
    object <- addCellMetadata(object, as.data.frame(enrichMatrix))

    if (plot){
    plotFeatures(object, target='cell', colorBy=colnames(enrichMatrix), ...)
    }
  
  return(object)
}

MySignaturesHeatmap <- function (object, topics = "all", selected.signatures = "all", 
    nCores = 4, aucMaxRank = 0.03 * nrow(aucellRankings), col.low = "dodgerblue", 
    col.mid = "floralwhite", col.high = "brown1", scale = TRUE, 
    ...) {
    # There was a typo at the time in the cisTopic source code that was preventing the usage of signatureHeatmap function
    # This is just the typo-corrected version

    if ((length(object@signatures) < 1)&&(selected.signatures!='annotation')) {
        stop("Please, run getSignaturesRegions() first.")
    }
    if (!"fastcluster" %in% installed.packages()) {
        stop("Please, install fastcluster: \n install.packages(\"fastcluster\")")
    }
    else {
        require(fastcluster)
    }
    if (!"ComplexHeatmap" %in% installed.packages()) {
        stop("Please, install ComplexHeatmap: source(\"https://bioconductor.org/biocLite.R\") \nbiocLite(\"ComplexHeatmap\")")
    }
    else {
        require(ComplexHeatmap)
    }
    scores <- .getScores(object)
    if (selected.signatures[1] != "annotation") {
        signatures <- object@signatures
        if (is.null(signatures)) {
            stop("Please run getSignaturesRegions() first.")
        }
        if (selected.signatures[1] != "all") {
            if (sum(selected.signatures %in% names(signatures)) != 
                length(selected.signatures)) {
                stop("Check whether the selected signatures have been stored in object@signatures.")
            }
            signatures <- signatures[selected.signatures]
        }
    }
    else {
        signatures <- split(object@region.data, object@region.data$annotation)
        if (is.null(signatures)) {
            stop("Please run annotateRegions() first.")
        }
        signatures <- lapply(signatures, function(x) rownames(x))
    }
    if (topics[1] != "all") {
        scores <- scores[, topics]
    }
    aucellRankings <- AUCell_buildRankings(as.matrix(scores), 
        nCores = nCores, plotStats = FALSE, verbose = FALSE)
    modulesAUC <- AUCell_calcAUC(signatures, aucellRankings, 
        nCores = nCores, aucMaxRank = aucMaxRank, verbose = FALSE)
    enrichMatrix <- getAUC(modulesAUC)
    if (scale) {
        enrichMatrix <- t(scale(t(enrichMatrix)))
        name_heatmap <- "Normalised AUC score"
    }
    else {
        name_heatmap <- "AUC score"
    }
    cl.topics <- fastcluster::hclust.vector(t(enrichMatrix), 
        method = "ward", metric = "euclidean")
    dd.col <- as.dendrogram(cl.topics)
    colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, 
        col.high))
    heatmap <- ComplexHeatmap::Heatmap(data.matrix(enrichMatrix), 
        col = colorPal(20), cluster_columns = dd.col, show_column_names = TRUE, 
        show_row_names = TRUE, heatmap_legend_param = list(legend_direction = "horizontal", 
            legend_width = unit(5, "cm"), title_position = "topcenter"), 
        name = name_heatmap, column_title_gp = gpar(fontface = "bold"), 
        ...)
    ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom")
}

region_annotate <- function(model, output_path=NULL, plot=NULL){
    # Annotate the cells with a score genomic regions (promoters, enhancers, ...) based on the cell-region distribution

    if(!'TxDb.Mmusculus.UCSC.mm10.knownGene' %in% installed.packages()){
        stop('Please install UCSV packages first')
    } else{
        require('TxDb.Mmusculus.UCSC.mm10.knownGene')
        txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
    }
    if(!is.null(plot)){
        pdf(paste(plot, model@project.name, '_region_annotation.pdf', sep=''))
    }
    model <- annotateRegions(model, txdb = txdb, annoDb = 'org.Mm.eg.db')
    model <- getRegionsScores(model, method = 'NormTop', scale=TRUE)
    model <- binarizecisTopics(model, thrP=0.975, plot=F)
    pred.matrix <- predictiveDistribution(model)
    aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)
    model <- MySignatureCellEnrichment(model, aucellRankings, selected.signatures='annotation', aucMaxRank = 0.1*nrow(aucellRankings), plot=FALSE)
    colnames(model@cell.data) <- make.names(colnames(model@cell.data), unique=TRUE)
    signaturesHeatmap(model, selected.signatures = 'annotation')
    if(!is.null(output_path)){
            saveRDS(model, file=paste(output_path, model@project.name, '_region_annotated.Rds', sep=''))
    }
    if(!is.null(plot)){
        dev.off()
    }
    return(model)
}

.getScores <- function (object) {
    scores <- object@region.data[, grep("Scores_Topic", colnames(object@region.data))]
    if (ncol(scores) < 1) {
        stop("Please, run getRegionsScores() first.")
    }
    colnames(scores) <- paste("Topic", 1:ncol(scores), sep = "")
    return(scores)
}

motifScoring <- function(model, path_to_bed, save_path){
    # Subset the scATAC binary accessibility matrix from a CisTopic model based on bed files (motif bed files)
    # Produces a matrix of size (number of motifs)x(number of cells)
    
    if (!'data.table' %in% installed.packages()){
        stop('Please install data.table')
    }
    else{
        require(data.table)
    }       
    
    binary_matrix <- model@count.matrix
    regions <- strsplit(row.names(binary_matrix),"[:-]")
    regions <- transpose(as.data.table(regions));
    names(regions) <- c('chr', 'start', 'end')
    # rownames(regions) <- c()
    regions$start <- strtoi(regions$start)
    regions$end <- strtoi(regions$end)
    
    # Creating the Granges object with the binary accessibility as metadata
    print('Creating GRanges object from the scATAC data')
    scatac_range <- GRanges(seqnames = Rle(regions$chr),
                      ranges = IRanges(start = regions$start, end = regions$end),
                      mat = as.matrix(binary_matrix))
    
    print('Done')
    # Looping through the bed files
    file.names <- dir(path_to_bed, pattern='.bed')
    
    #Initialize the motif scoring dataframe
    motif_scores <- data.frame(matrix(nrow=length(file.names),
                                      ncol=dim(binary_matrix)[2]))
    print('Initialise')
    print(dim(motif_scores))
    for (i in 1:length(file.names)){
        print(i)
        motif_range <- import(paste(path_to_bed, file.names[i], sep=''))
        sm <- subsetByOverlaps(x = scatac_range, ranges = motif_range)
        scores <- colSums(as.matrix(sm@elementMetadata))
        motif_scores[i,] <- as.data.frame(t(scores))
        rownames(motif_scores)[i] <- gsub(pattern='.bed', '', file.names[i])
    }
    colnames(motif_scores) <- colnames(model@count.matrix)
    saveRDS(motif_scores, file=paste(save_path, 'motif_scores.Rda', sep=''))
    print('Done!')
}

find_k_neighbors = function(k=20, dist.matrix){
  # Find the k closest neighbors from cells using a distance matrix 

  require(FastKNN)

  n = dim(dist.matrix)[1] # number of WT cells 
  nn = matrix(0, n, k)
  
  for (i in 1:n){
    nn[i, ] = k.nearest.neighbors(i, dist.matrix, k=k)
  }
  
  .get_col_name <- function(index){return(colnames(dist.matrix)[index])}
  
  nn_names = apply(nn, MARGIN = c(1,2), FUN=.get_col_name)
  rownames(nn_names) = rownames(dist.matrix)
  return(nn_names)
}

.get_TF_score = function(cell_name, df_table, transcription_factor){
    #helper function of average_TF_score
    if (transcription_factor%in%rownames(df_table)){
        if(cell_name%in%colnames(df_table)){
            return(df_table[transcription_factor, cell_name])
    } else {
        print(paste(cell_name, transcription_factor))
        stop('Invalid tables')}
    } 
}

average_TF_score = function(df.neighbors, deviation_table){
    # Compute the average of TF deviation scores over the neighbors of reference cells
    # Inputs : 
    # -> df.neighbors (reference cells X neighbors)
    # -> deviation_table : table of TF scores from chromVAR for all the neighbors of reference cells (TF X cells)

    df = data.frame(matrix(nrow=dim(df.neighbors)[1], ncol=dim(deviation_table)[1]))
    colnames(df) = rownames(deviation_table)
    rownames(df) = rownames(df.neighbors)
    for (i in 1:nrow(df.neighbors)){
        vec = rowMeans(deviation_table[, df.neighbors[i, ], drop=FALSE])
        df[i, ] = vec
    }
    return(df)
}

quantile_normalization = function(df){
  df_rank = apply(df,2,rank,ties.method='min')
  df_sorted = data.frame(apply(df, 2, sort))
  df_mean = apply(df_sorted, 1, mean)
  
  index_to_mean = function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Function for downsampling
downsamp_one = function(x,n){
  hist(sample(rep(1:length(x),x),size = n,replace=F),breaks=(0:length(x)+.5),plot=F)$counts
}
downsample = function(u,min_umis){
  cell_mask=colnames(u)[Matrix::colSums(u,na.rm=T)>min_umis]  
  message("Downsampling ", length(cell_mask), " cells to ",min_umis," UMIs")
  ds=as.matrix(apply(u[,cell_mask],2,downsamp_one,min_umis))
  rownames(ds)=rownames(u)
  return(ds)
}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right"), top='') {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1, top=top,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2, top=top,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined) 
}

getCountsFromFrags <- function(fragment_file_XX, peaks_gr, barcodes){
    # Adapted code from Caleb Lareau's github https://github.com/caleblareau/scATAC_10X_raw_to_kmers/blob/master/example_kmers.R
    #
    # Create SummarizedExperiment object from 10X scATAC-seq fragment file 
    # given a set of peaks (peaks_gr) and filtered barcodes (barcodes)

    require(GenomicRanges)
    require(SummarizedExperiment)
    require(data.table)
    require(dplyr)

    # GRanges of all the reads that are in the filtered cells 
    fragments = data.table::fread(fragment_file_XX) %>%
        data.frame() %>% filter(V4 %in% barcodes) %>%
        GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)

    # Get nb of reads 
    denom <- table(GenomicRanges::mcols(fragments)$V4)
    barcodes_found <- names(denom)

    ovPEAK <- GenomicRanges::findOverlaps(peaks_gr, fragments)

    id <- factor(as.character(GenomicRanges::mcols(fragments)$V4), levels = barcodes_found)

    # Make sparse matrix with counts with peaks by  unique barcode
    countdf <- data.frame(peaks = S4Vectors::queryHits(ovPEAK),
                        sample = as.numeric(id)[S4Vectors::subjectHits(ovPEAK)]) %>%
    dplyr::group_by(peaks,sample) %>% dplyr::summarise(count = n()) %>% data.matrix()
  
    m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks_gr)),
                            j = c(countdf[,2], length(barcodes_found)),
                            x = c(countdf[,3],0))
    
    colnames(m) <- barcodes_found
      
    # Make a polished colData
    colData <- data.frame(
    sample = barcodes_found,
    depth = as.numeric(denom),
    FRIP = Matrix::colSums(m)/as.numeric(denom)
    )
    # Make sure that the SE can be correctly constructed
    stopifnot(all(colData$sample == colnames(m)))

    # Make summarized Experiment
    SE <- SummarizedExperiment::SummarizedExperiment(
    rowRanges = peaks_gr,
    assays = list(counts = m),
    colData = colData
    )
  return(SE)
}

format_bed_to_homer <- function(bedfile){
    # format a bed file witg 3 columns to HOMER format

    require(readr)

    print('Reading file...')
    table = read.delim(bedfile, header=FALSE)
    colnames(table) = c('chr', 'start', 'end')
    table$PeakID = 1:nrow(table)
    table$Unused = ''
    table$strand = '.'

    path_to_save = gsub(bedfile, pattern='.bed', replacement='_HOMER_format.bed')
    print(path_to_save)

    write.table(table, file=path_to_save, quote=FALSE,col.names=FALSE, row.names=FALSE, sep='\t')
}
               
createcistopic = function(path_to_XX, sample_list, min_cells, name, save_path){
    # Initialize cistopic object and filter cells and peaks

    require(cisTopic)
    require(Matrix)

    matrix = Matrix::readMM(paste0(path_to_XX, 'matrix.mtx'))
    peaks = read.table(paste0(path_to_XX, 'peaks.bed'), header=FALSE)
    barcodes = read.table(paste0(path_to_XX, 'barcodes.tsv'), header=FALSE)

    colnames(matrix) = barcodes$V1
    peaks$site_name = paste0(peaks$V1, ':', peaks$V2, '-', peaks$V3)
    rownames(matrix) = peaks$site_name

    matrix_cells = matrix[, sample_list]

    bmat = matrix_cells
    bmat@x[bmat@x > 0] = 1 # binarize matrix to count coverage per peak

    coverage = Matrix::rowSums(bmat)
    indices = coverage > min_cells

    count.matrix = matrix_cells[indices, ]

    cistopic = createcisTopicObject(count.matrix=count.matrix,
                                    project.name=name)

    saveRDS(cistopic, save_path)

    return(cistopic)
}
