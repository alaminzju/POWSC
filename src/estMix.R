######################################################################
######### Using Seurat to estimate the mixture probabilities ######### 
######################################################################

seuratMix = function(counts){
    require(Seurat)
    seuset = CreateSeuratObject(raw.data = counts,
        min.cells = 3, min.genes = 200)
    seuset = NormalizeData(object = seuset, 
        normalization.method = "LogNormalize", 
        scale.factor = 10000)
    seuset = FindVariableGenes(object = seuset,
        mean.function = ExpMean, dispersion.function = LogVMR, 
        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    seuset = ScaleData(object = seuset)
    seuset = RunPCA(object = seuset, 
        pc.genes = seuset@var.genes, o.print = TRUE, 
        pcs.print = 1:5, genes.print = 5)
    seuset = FindClusters(object = seuset, 
        reduction.type = "pca", dims.use = 1:8, 
        resolution = 1.0, print.output = 0, 
        save.SNN = TRUE)
    return(seuset)
}


#######################################################################
########### Using SC3 to estimate the mixture probabilities ########### 
#######################################################################

sc3Mix = function(counts, ks = NULL){
    require(SC3)
    sce = SingleCellExperiment(
        assays = list(
            counts = as.matrix(counts),
            logcounts = log2(as.matrix(counts) + 1)
        )
    )
    rowData(sce)$feature_symbol <- rownames(sce)
    sce = sce[!duplicated(rowData(sce)$feature_symbol), ]
    sce = sc3_prepare(sce)
    if(is.null(ks)){
        sce = sc3_estimate_k(sce)
        cat("estimated k clusters: ", metadata(sce)$sc3$k_estimation, "\n")
        k = metadata(sce)$sc3$k_estimation
        ks = c(k-1, k, k+1)
    }else{
        ks = ks
    }
    sce = sc3_calc_dists(sce)
    sce = sc3_calc_transfs(sce)
    sce = sc3_kmeans(sce, ks = ks)
    sce = sc3_calc_consens(sce)
    colTb = data.frame(colData(sce))
    nk = length(ks)
    pList = lapply(seq_len(nk), function(i){
        tb = table(colTb[,i])/ sum(table(colTb[,i]))
        tb = sort(tb, decreasing = T)
        return(tb)
    })
    names(pList) = paste0(ks, "_clusters ")
    return(list(table = colTb, percentage = pList))
}






