########################################################
## Power function, compute the stratified power-realted quantities
##
## Return values:
## TD: number of True Discoveries in each stratum
## FD: number of False Discoveries in each stratum
## power: within strata, proportion of TD out of total DE
## alpha.nomial: cutoff of alpha on raw p-values
## alpha: empirical p-value in each stratum
## alpha.marginal: overall empirical p-value
########################################################
PowerEst = function(fdr, alpha, Zg, Zg2, xgr){
    ## p is input nominal p-value or FDR.
    ## alpha is cutoff of p
    ## !Zg is indicator for TN,  Zg2 is indicators for TP,
    ## xgr is grouping in covariate

    ix.D = fdr <= alpha
    N = sum(ix.D) ## this is the total discovery
    N.stratified = tapply(ix.D, xgr, sum)

    ##  TD
    id.TP = Zg2==1
    TD = tapply(fdr[id.TP] <= alpha, xgr[id.TP], sum)
    TD[is.na(TD)] = 0

    ##  FD
    id.TN = Zg==0
    FD = tapply(fdr[id.TN] <= alpha, xgr[id.TN], sum)
    FD[is.na(FD)] = 0

    ## type I error
    alpha = as.vector(FD/table(xgr[id.TN]))
    alpha.marginal = sum(FD)/sum(id.TN)

    ## power & precision
    power=as.vector(TD/table(xgr[id.TP]))
    power.marginal=sum(TD,na.rm=TRUE)/sum(id.TP)

    ## FDR
    FDR = FD / N.stratified
    FDR.marginal = sum(FD, na.rm=TRUE) / N
    
    ## conditional truths
    CT = table(xgr[id.TP])
    
    list(TD=TD,FD=FD,alpha.nominal=alpha, CT = CT,
         alpha=alpha, alpha.marginal=alpha.marginal,
         power=power, power.marginal=power.marginal,
         FDR=FDR, FDR.marginal=FDR.marginal)
}


## Continous case corresponding to the Phase II DE, delta means lfc
Power_Cont = function(DErslt, simData, alpha = 0.1, delta = 0.5, strata = c(0,10,2^(1:4)*10,Inf)){
    fdrvec = DErslt$cont$fdr
    lfc = simData$lfc
    ngenes = nrow(simData$sce)
    DEid = simData$ix.DE2
    Zg = Zg2 = rep(0, ngenes)
    Zg[DEid] = 1
    ix = which(abs(lfc) > delta)
    Zg2[DEid[ix]] = 1
    sce = simData$sce
    Y = round(assays(sce)[[1]])
    # sizeF = colSums(Y)
    # sizeF = sizeF/median(sizeF)  
    # X.bar = rowMeans(sweep(Y,2,sizeF,FUN="/"))
    X.bar = rowMeans(Y)
    ix.keep = which(X.bar>0)
    xgr = cut(X.bar[ix.keep], strata)
    
    # lev = levels(xgr)
    # ix.keep = ix.keep[!(xgr %in% lev[strata.filtered])]
    # ## recut
    # xgr = cut(X.bar[ix.keep], strata[-strata.filtered])
    
    # Interested genes
    Zg = Zg[ix.keep]; Zg2 = Zg2[ix.keep]
    fdrvec = fdrvec[ix.keep]
    
    pow = PowerEst(fdrvec, alpha, Zg, Zg2, xgr=xgr)
    return(pow)
}


## Discreate case corresponding to the Phase I DE, delta means pi.df
Power_Disc = function(DErslt, simData, alpha = 0.1, delta = 0.1, strata = seq(0, 1, by = 0.2)){
    fdrvec = DErslt$disc$fdr
    pi.df = simData$pi.df
    ngenes = nrow(simData$sce)
    DEid = simData$ix.DE1
    Zg = Zg2 = rep(0, ngenes)
    Zg[DEid] = 1
    ix = which(abs(pi.df) > delta)
    Zg2[DEid[ix]] = 1
    sce = simData$sce
    ntotal = ncol(sce)
    Y = round(assays(sce)[[1]])
    no0rate = 1 - apply(Y, 1, function(x) sum(x ==0))/ntotal
    ix.keep = which(no0rate > 0.01) # too small none 0 rate cannot detect
    xgr = cut(no0rate[ix.keep], strata)
    
    # lev = levels(xgr)
    # ix.keep = ix.keep[!(xgr %in% lev[strata.filtered])]
    # ## recut
    # xgr = cut(X.bar[ix.keep], strata[-strata.filtered])
    
    # Interested genes
    Zg = Zg[ix.keep]; Zg2 = Zg2[ix.keep]
    fdrvec = fdrvec[ix.keep]
    
    pow = PowerEst(fdrvec, alpha, Zg, Zg2, xgr=xgr)
    return(pow)
}






