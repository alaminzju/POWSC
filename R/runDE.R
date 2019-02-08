#### RUN DE METHODS #####
#@param rawEset much has a pData variable called celltype indicator

library(Biobase)
library(SC2P)
# Set up global variables
model = ~ cellTypes + cngeneson # in MAST
termTotest = "cellTypes" # in SC2P


toEset = function(sce){
    y = assays(sce)$counts
    cellTypes = colData(sce)$cellTypes
    phenoData = new("AnnotatedDataFrame", data = data.frame(cellTypes = cellTypes))
    eset = ExpressionSet(assayData = y, phenoData = phenoData)
    return(eset)
}
## SC2P
runSC2P = function(sce) {
    rawEset = toEset(sce = sce); ngene = nrow(sce)
    norm = eset2Phase(rawEset)  ## eset if of class ExpressionSet or eSet
    deSC2P = twoPhaseDE(norm, design = termTotest, test.which=1, offset="sf")
    discPval = deSC2P$Ph1.pval; discFdr = p.adjust(discPval, method = "BH")
    contPval = deSC2P$Ph2.pval; contFdr = p.adjust(contPval, method = "BH")
    disc = data.frame(geneIndex = 1:length(discPval), pval = discPval, fdr = discFdr)
    cont = data.frame(geneIndex = 1:length(contPval), pval = contPval, fdr = contFdr)
    rownames(disc) = rownames(cont) = rownames(deSC2P)
    disc = disc[complete.cases(disc), ]; cont = cont[complete.cases(cont), ]
    return(list(table = deSC2P, cont = cont, disc = disc))
}

## MAST
runMAST = function(sce) {
    ## grab raw counts and compute log TPM
    rawEset = toEset(sce = sce)
    Y = exprs(rawEset)
    sf = colSums(Y)/1e6
    ltpm = log2(sweep(Y, 2, sf, FUN="/")+1)
    sca = FromMatrix(ltpm, pData(rawEset), fData(rawEset)) # the object is sca for MAST
    ## Start testing
    cdr2 = colSums(assay(sca)>0)
    colData(sca)$cngeneson = scale(cdr2)
    thres = thresholdSCRNACountMatrix(assay(sca)) #, nbins=200, min_per_bin=30)
    assays(sca) = list(thresh=thres$counts_threshold, tpm=assay(sca))
    fit = zlm(model, sca)
    rslt = lrTest(fit, termTotest)  # fiting zlm model
    table = rslt[,,3]
    two.des = list()
    for (i in 1:2){
        pval = table[, i]
        fdr = p.adjust(pval, method = "BH")
        result = data.frame(geneIndex=1:length(pval), pval=pval, fdr=fdr)
        res = result[complete.cases(result),]
        two.des[[i]] = res
    }
    return(list(table = table, cont = two.des[[1]], disc = two.des[[2]]))
}


## runDE summary function 
runDE = function(sce, DE_Method = c("MAST", "SC2P")){
    DE_Method = match.arg(DE_Method)
    if (DE_Method == "MAST"){
        DE_rslt = runMAST(sce)
    }else{
        DE_rslt = runSC2P(sce)
    }
    return(DE_rslt)
}





##############################################################################
############################## other DE methods ############################## 
##############################################################################
## BPSC
runBPSC = function(rawEset) {
    require(BPSC)
    Y = exprs(rawEset)
    sf = colSums(Y)/1e6
    bp_mat = sweep(Y, 2, sf, FUN="/")+1
    groups = as.factor(droplevels(pData(rawEset)$celltype))
    controlIds = which(groups==levels(groups)[1])
    design = model.matrix( ~ groups)
    rslt = BPglm(data = bp_mat,
                     controlIds = controlIds, design = design, coef = 2,
                     estIntPar = FALSE)
    table = summary(rslt)$topTable
    pval = table[, 3]
    pval[is.na(pval)] <- 1
    qval  = qvalue(pval)$qvalues

    ## construct results, ordered by the orignal order of genes
    ix = sort(pval, index.return=TRUE)$ix
    result = data.frame(geneIndex=ix, pval=pval[ix], qval=qval[ix])
    res = result[order(result$geneIndex),]
    return(res)
}

# SCDD
runSCDD = function(simData){
    require("SingleCellExperiment"); require("scDD")
    Y = exprs(rawEset)
    sf = colSums(Y)/1e6
    ltpm = log2(sweep(Y, 2, sf, FUN="/")+1)
    group = as.factor(droplevels(pData(rawEset)$celltype))
    priorParam = list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
    sce = SingleCellExperiment(assays=list(normcounts = as.matrix(ltpm)), colData=data.frame(condition = group))
    rslt = scDD(sce, priorParam, testZeroes=F)
    table = results(rslt)
    pval = table[, 7]
    pval[is.na(pval)] <- 1
    fdr  = qvalue::lfdr(pval)

    ## construct results, ordered by the orignal order of genes
    ix = sort(pval, index.return=TRUE)$ix
    result = data.frame(geneIndex=ix, pval=pval[ix], fdr=fdr[ix])
    res = result[order(result$geneIndex),]
    return(res)
}


## SCDE
runSCDE = function(rawEset, ncores=32) {
    require(scde)
    Y = exprs(rawEset)
    Y = apply(Y,2,function(x) {
        storage.mode(x) <- 'integer'; x})
    colnames(Y) = 1:ncol(Y)
    cd = clean.counts(Y)
    sg = (pData(rawEset)$celltype)[as.integer(colnames(cd))]
    o.ifm = scde.error.models(counts=cd, groups=sg, n.cores=ncores,
                               threshold.segmentation = TRUE,
                               save.crossfit.plots = FALSE,
                               save.model.plots = FALSE, verbose = 1)

    valid.cells = o.ifm$corr.a > 0
    table(valid.cells)
    o.ifm = o.ifm[valid.cells, ]
    o.prior = scde.expression.prior(models = o.ifm,
                                     counts = cd, length.out = 400,
                                     show.plot = FALSE)
    groups = sg[valid.cells]
    de.scde = scde.expression.difference(o.ifm, cd, o.prior,
                                          groups  =  groups,
                                          n.randomizations  =  100,
                                          n.cores  =  ncores,
                                          verbose  =  0)
    return(de.scde)
}


makeROC = function(pval, trueDE) {
    pred = prediction(-log10(pval), trueDE)
    auc = performance(pred,"auc")
    auc = unlist(slot(auc, "y.values"))
    perf = performance(pred,"tpr","fpr")
    return(list(perf = perf, auc = auc))
}


### Main function #### should rewrite this function
# runDE = function(simData){
#     rawEset = simData$eset
#     rslt = runMAST(rawEset)
#     pvals = rslt$pval
#     trueDE = sign(match(rownames(rslt), simData$DEGs, nomatch = 0))
#     roc_rslt = makeROC(rslt$pval, trueDE)
#     return(list(simData = simData,
#                 auc = roc_rslt$auc, perf = roc_rslt$perf))
# }


