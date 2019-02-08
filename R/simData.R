rzip <- function(n, p0, lambda) {
    y <- rep(0, n)
    ix <- rbinom(n, size=1, prob=p0)
    y[!ix] <- rpois(sum(!ix), lambda)
    y
}

### generate LNP random variable
## The real data is not exactly LNP, the right tail is much shorter
## Using LNP to simulate there are too many large numbers.
## How to shave the right tail??
rLNP <- function(n, mu, sigma, sf) {
    theta <- 2^rnorm(n, mean=mu, sd=sigma)
    ## should shave the right tail of theta a bit to avoid excessive large number??

    y <- suppressWarnings(rpois(n, theta*sf))
    if (sum(is.na(y)) > 0){
        m = max(y[!is.na(y)])
        y[is.na(y)] = m
    }
    return(y)
}

###############################################
################ Main Function ################
###############################################
GenerateCountMatrix <- function(pi.g, p0, lambda, mu, sigma, sf){

    stopifnot(identical(length(p0), length(lambda), length(sf)),
              identical(length(pi.g), length(mu), length(sigma)))

    N <- length(p0)
    G <- length(mu)

    ## simulate Z, indicator for being in foreground.
    Z <- matrix(rbinom(N*G, size=1, prob=pi.g), ncol=N)

    ## loop over cells to generate counts
    Y <- matrix(0, nrow=G, ncol=N)
    for(i in 1:N) {
        ix <- Z[,i] == 0
        Y[ix, i] <- rzip(sum(ix), p0[i], lambda[i])
        Y[!ix, i] <- rLNP(sum(!ix), mu[!ix], sigma[!ix], sf[i])
    }

    return(Y)
}

######################################################################
########## Function 1 use one cell type to simulate another ########## 
######################################################################
# SimulateEset = function(n1 = 100, n2 = 100, perDE = 0.05, estParas) {
#     # Parameters from estimations
#     pi.g = estParas$pi.g
#     p0 = estParas$p0
#     lambda = estParas$lambda
#     mu = estParas$mu
#     sigma = estParas$sd
#     sf = estParas$sf
# 
#     ngene = length(pi.g)
#     nDE1 = nDE2 = ngene * perDE
#     ## simulate index for DE genes.
#     ## I'm picking from genes with large means.
#     ## If DE genes are for genes with very small mean, they won't be detected.
#     n0 = sum(mu >= 3)
#     ix.highGenes = order(mu, decreasing=TRUE)[1:n0]
#     ix.DE1 = sample(ix.highGenes, nDE1)
#     ix.DE2 = sample(ix.highGenes, nDE2)
# 
#     ix.DEs = union(ix.DE1, ix.DE2)
#     ## general parameters in both groups
#     sf1 = sample(sf, n1, replace=TRUE)
#     sf2 = sample(sf, n2, replace=TRUE)
#     p0.1 = sample(p0, n1, replace=TRUE )
#     p0.2 = sample(p0, n2, replace=TRUE)
#     lambda1 = sample(lambda, n1, replace=TRUE)
#     lambda2= sample(lambda, n2, replace=TRUE)
# 
#     ## generate parameters for second group
#     ## the inputs are assumed to be for first group
#     ## for mixing proportion
#     
#     ## for zero ratio in phase I 
#     pi.g2 = pi.g
#     tmp = pi.g[ix.DE1]
#     tmp[tmp<0.5] = tmp[tmp<0.5] + runif(sum(tmp<0.5), 0.1, 0.3)
#     tmp[tmp>=0.5] = tmp[tmp>=0.5] - runif(sum(tmp>=0.5), 0.1, 0.3)
#     pi.g2[ix.DE1] = tmp
#     pi.df = pi.g2 - pi.g
#     
#     ## for mean expression in phase II
#     mu2 = mu
#     tmp = c(rnorm(1000, mean=-1, sd=1), rnorm(1000, mean=1, sd=1))
#     mu.diff = sample(tmp, length(ix.DE2))
#     mu2[ix.DE2] = mu[ix.DE2] + mu.diff
#     ## I will not change the variance parameters (sigma or phi will be the same in two groups)
# 
#     ## For Phase II DEGs lfc
#     lfc  = mu.diff
#     ## generate counts in two groups
#     y1 = GenerateCountMatrix(pi.g, p0.1, lambda1, mu, sigma, sf1)
#     y2 = GenerateCountMatrix(pi.g2, p0.2, lambda2, mu2, sigma, sf2)
#     y = cbind(y1, y2)
#     rownames(y) = paste0("g", 1:nrow(y))
#     celltype = rep(paste0("celltype", c(1,2)), c(n1, n2))
#     phenoData = new("AnnotatedDataFrame", data = data.frame(celltype = celltype))
#     eset = ExpressionSet(assayData = y,
#                           phenoData = phenoData)
#     DEGs = paste0("g", union(ix.DE1, ix.DE2))
#     list(ix.DE1=ix.DE1, ix.DE2=ix.DE2, ix.DEs = ix.DEs, DEGs = DEGs, eset=eset,
#          pi.g1=pi.g, pi.g2=pi.g2, mu1=mu, mu2=mu2, lfc = lfc, pi.df=pi.df, ngenes = ngene)
# }
# 






##################################################################
####### Function 2 considering two (pair-wised) conditions ####### 
##################################################################
Simulate2SCE = function(n = 100, perDE = 0.05, estParas1, estParas2) {
    # equally divide the sample size 
    n1 = n2 = round(n / 2)
    # Parameters from estimations
    pi.g1 = estParas1$pi.g; pi.g2 = estParas2$pi.g
    p01 = estParas1$p0; p02 = estParas2$p0
    lambda1 = estParas1$lambda; lambda2 = estParas2$lambda
    mu1 = estParas1$mu; mu2 = estParas2$mu
    sigma1 = estParas1$sd; sigma2 = estParas2$sd
    sf1 = estParas1$sf; sf2 = estParas2$sf
    
    if (length(pi.g1) == length(pi.g2)){
        ngene = length(pi.g1)
    }else{
        stop("Please Input Proper Parameter Estimation Objects")
    }
    ## simulate index for DE genes.
    ## I'm picking from genes with large means.
    ## If DE genes are for genes with very small mean, they won't be detected.
    n0 = max(sum(mu1 > 3), sum(mu2 > 3))
    nDE1 = nDE2 = ngene * perDE
    ix1.highGenes = order(mu1, decreasing=TRUE)[1:n0]
    ix2.highGenes = order(mu2, decreasing=TRUE)[1:n0]
    
    ix.DE1 = sample(union(ix1.highGenes, ix2.highGenes), nDE1)
    ix.DE2 = sample(union(ix1.highGenes, ix2.highGenes), nDE2)
    
    ix.DEs = union(ix.DE1, ix.DE2)
    ## general parameters in both groups
    sf1 = sample(sf1, n1, replace=TRUE)
    sf2 = sample(sf2, n2, replace=TRUE)
    p0.1 = sample(p01, n1, replace=TRUE )
    p0.2 = sample(p02, n2, replace=TRUE)
    lambda1 = sample(lambda1, n1, replace=TRUE)
    lambda2= sample(lambda2, n2, replace=TRUE)
    
    ## For Phase I DEGs zero ratio
    tmp = pi.g2[ix.DE1]  
    tmp[tmp < 0.5] = tmp[tmp < 0.5] + runif(sum(tmp<0.5), 0.1, 0.3)
    tmp[tmp >= 0.5] = tmp[tmp >= 0.5] - runif(sum(tmp >=0.5), 0.1, 0.3)
    pi.g2[ix.DE1] = tmp
    pi.df = tmp - pi.g1[ix.DE1]
    
    ## For Phase II DEGs lfc
    tmp = c(rnorm(1000, mean=-1, sd=1), rnorm(1000, mean=1, sd=1))
    mu.diff = sample(tmp, length(ix.DE2))
    mu2[ix.DE2] = mu1[ix.DE2] + mu.diff
    lfc  = mu.diff
    
    ## generate counts in two groups
    y1 = GenerateCountMatrix(pi.g1, p0.1, lambda1, mu1, sigma1, sf1)
    y2 = GenerateCountMatrix(pi.g2, p0.2, lambda2, mu2, sigma2, sf2)
    y = cbind(y1, y2)
    rownames(y) = paste0("g", 1:nrow(y))
    celltypes = rep(paste0("celltype", c(1,2)), c(n1, n2))
    sce = SingleCellExperiment(
        assays = list(counts = y),
        rowData = data.frame(geneNames = rownames(y), stringsAsFactors = F),
        colData = data.frame(cellTypes = celltypes, stringsAsFactors = F)
    )
    DEGs = paste0("g", union(ix.DE1, ix.DE2))
    list(ix.DE1=ix.DE1, ix.DE2=ix.DE2, ix.DEs = ix.DEs, DEGs = DEGs, sce = sce,
         pi.g1=pi.g1, pi.g2=pi.g2, mu1=mu1, mu2=mu2, lfc = lfc, pi.df=pi.df, ngenes = ngene)
}






##################################################################
##### Function 3 considering a mixture cell type conditions ######
##################################################################
SimulateMultiSCEs = function(n = 100, perDE = 0.05, estParas, multiProb) {
    ## Parameters from estimations
    pi.g = estParas$pi.g
    p0 = estParas$p0
    lambda = estParas$lambda
    mu = estParas$mu
    sigma = estParas$sd
    sf = estParas$sf
    ngene = length(pi.g)
    ncelltype = length(multiProb)
    ## Simulate index for DE genes. Picking from genes with large means otherwise 
    ## genes with very small mean won't be detected.
    nhighGenes = sum(mu > 3)
    ix.highGenes = order(mu, decreasing=TRUE)[1:nhighGenes]
    n0 = nhighGenes * perDE
    
    ## Simulate number of the cell types from multinomial distribution
    
    multiSize = as.vector(rmultinom(1, n, multiProb))
    simAll = lapply(multiSize, function(tmpSize){
        ix.DE1 = sample(ix.highGenes, n0)
        ix.DE2 = sample(ix.highGenes, n0)
        ix.DEs = union(ix.DE1, ix.DE2)
        ## general new parameters in one new cell type
        tmp.sf = sample(sf, tmpSize, replace=TRUE)
        tmp.p0 = sample(p0, tmpSize, replace=TRUE )
        tmp.lambda = sample(lambda, tmpSize, replace=TRUE)
        ## perturbation of pi.g / zero ratio in phase I 
        tmp.pi.g = pi.g
        tmp = pi.g[ix.DE1]
        tmp[tmp<0.5] = tmp[tmp<0.5] + runif(sum(tmp<0.5), 0.1, 0.3)
        tmp[tmp>=0.5] = tmp[tmp>=0.5] - runif(sum(tmp>=0.5), 0.1, 0.3)
        tmp.pi.g[ix.DE1] = tmp
        pi.df = tmp.pi.g - pi.g
        ## perturbation of mu / lfc in phase II
        tmp.mu = mu
        tmp = c(rnorm(1000, mean=-1, sd=1), rnorm(1000, mean=1, sd=1))
        mu.diff = sample(tmp, length(ix.DE2))
        tmp.mu[ix.DE2] = mu[ix.DE2] + mu.diff
        lfc  = mu.diff
        ## generate counts in two groups
        tmp.y = GenerateCountMatrix(tmp.pi.g, tmp.p0, tmp.lambda, tmp.mu, sigma, tmp.sf)
        tmp.rslt = list(Y = tmp.y, ix.DE1 = ix.DE1, ix.DE2 = ix.DE2, pi.g = tmp.pi.g, mu = tmp.mu)
        return(tmp.rslt)
    })
    
    SimulateMultiSCEs = NULL
    for (i in 1:ncelltype) {
        for (j in i: ncelltype) {
            if (j  > i){
                comID = c(i, j)
                comName = paste0(i, "_vs_", j)
                Y = do.call(cbind, lapply(simAll[comID], function(x) x$Y))
                ix.DE1 = do.call(union, lapply(simAll[comID], function(x) x$ix.DE1))
                ix.DE2 = do.call(union, lapply(simAll[comID], function(x) x$ix.DE2))
                pi.g = do.call(cbind, lapply(simAll[comID], function(x) x$pi.g))
                pi.df = (pi.g[,1] - pi.g[,2])[ix.DE1]
                mu = do.call(cbind, lapply(simAll[comID], function(x) x$mu))
                lfc = (mu[,1] - mu[,2])[ix.DE2]
                rownames(Y) = paste0("g", 1:nrow(Y))
                celltypes = rep(paste0("celltype", comID), multiSize[comID])
                sce = SingleCellExperiment(
                    assays = list(counts = Y),
                    rowData = data.frame(geneNames = rownames(Y), stringsAsFactors = F),
                    colData = data.frame(cellTypes = celltypes, stringsAsFactors = F)
                )
                SimulateMultiSCEs[[comName]] = list(sce = sce, ix.DE1 = ix.DE1, ix.DE2 = ix.DE2, 
                                                pi.df = pi.df, lfc = lfc)
            }else{
                next
            }
        }
    }
    return(SimulateMultiSCEs)
}





