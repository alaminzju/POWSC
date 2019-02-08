## parameter estimation
RobustPoi0 <- function(x){
    tt=table(x)
    n0=sum(x==0)
    if (n0 == length(x)){
        c("p0"=1,"mu"=0)
    }else{
        if (names(tt)[1]=="0") {
            xx=as.numeric(names(tt))[-1]
            tt=as.vector(tt)[-1]
        }else{
            xx=as.numeric(names(tt))
            tt=as.vector(tt)
        }
        tt=log(tt)+log(gamma(1+xx))
        ##fit the regression without N0
        beta=lm(tt~xx,weight=1/exp(xx))$coef
        mu.hat=exp(beta[2])
        p0=(n0-exp(beta[1]))/(exp(beta[1]+mu.hat)+n0-exp(beta[1]))
        if (any(p0<0)) {warning("Proportion of zero inflation is negative")}
        p0 <- pmax(p0, 0)
        c("p0"=p0,"mu"=mu.hat)
    }
}


## Shrink the mu for sample sizeas are so small
shrink.mu=function(y,s,n){
    mu.g=rep(NA,length(y))
    k=which(n>1)
    if (length(k)<length(n)) {fill=TRUE} else {fill=FALSE}
    s=s[k];y=y[k];n=n[k]

    mu0=weighted.mean(y,w=n)

    s2.total=sum(s^2*(n-1))+sum(n*(y-mu0)^2)
    s2.total=s2.total/sum(n)

    s2.1=sum(s^2*(n-1))/sum(n)
    s2.0=s2.total-s2.1
    ### shrink mu
    mu.sub=  (y*n/s2.1+mu0/s2.0)/(n/s2.1+1/s2.0)
    mu.g[k]=mu.sub
    if (fill) mu.g[-k]=mu0

    mu.g
}

my.mad=function (x, center = median(x), constant = 1.4826, na.rm = FALSE){
    if (na.rm)
        x <- x[!is.na(x)]
    res=x-center
    constant * median(res[res>0])
}

dZinf.pois=function(x, p0, mu){
    (x==0)*(p0+(1-p0)*dpois(x,mu))+(x>0)*(1-p0)*dpois(x,mu)
}

dLNP2 <- function(x, mu, sigma, l=1){
    x.min <- pmax(0, x-0.5)
    pnorm(log2((x+0.5)/l), mu, sigma) - pnorm(log2(x.min/l), mu, sigma)
}



#########################################################
##################### Main Function #####################
#########################################################
Eset2Phase = function(sce, low.prob=0.99){ 
    # check the input whether it is an SingleCellExperiment Object or a matrix
    if (is(sce, "SingleCellExperiment")){
        Y = round(assays(sce)[[1]])
    }else{Y = round(sce)}
    ## ## initial estimate of prob(X\in bg)
    Cell0=colMeans(Y==0) # each cell has this percentage 0
    par1= apply(Y,2,function(yy) {
        yy=yy[yy<=15]
        RobustPoi0(yy)})
    pi0.hat=Cell0/(par1[1,]+(1-par1[1,])*dpois(0,par1[2,]))
    if (any((pi0.hat > 1))) {warning("Zero proportion is greater than estimation.")}
    pi0.hat <- pmin(pi0.hat, 1)
    prob0=pi0.hat*par1[1,]+ pi0.hat*(1-par1[1,])*dpois(0,par1[2,]) ## ZIP prob at 0
    ## First round
    ## get the 1-low.prob quantile of ZIP
    x0=qpois(pmax(1-(1-low.prob)/(1-par1[1,]),0),par1[2,])
    Z= sweep(Y,2,x0)>0 # indicate if a gene is > bg
    L=colSums(Y*Z)/1e6 # so far it is like simple total..

    mu.g1=log2(rowSums(Z*Y)/rowSums(sweep(Z,2,L,FUN="*")))
    mu.g1[is.na(mu.g1)]=0 ## if allZ is 0, it gets NA,
    ### but we should shrink mu.g1 as well since some mu.g1 is estimated by only a few observations
    ## leave it here for now.
    n.g1=rowSums(Z)
    y1=log2(sweep(Y,2,L,FUN="/")+1) #like TPM**
    s.g1=sqrt(rowSums(Z*sweep(y1,1,mu.g1)^2)/(n.g1-1)) ## TPM type of SD
    mu.g2 = shrink.mu(mu.g1,s.g1,n.g1)
    ## get sd.g
    res.g1=log2(sweep(Y,2,L,FUN="/")+1)-mu.g1
    ## mad of those res.g1 that are associated with Z==1
    tmp=array(0,dim=c(dim(res.g1),2))
    tmp[,,1]=res.g1;tmp[,,2]=Z
    sd.g1=apply(tmp,1,function(xx) my.mad(xx[xx[,2]==1,1]))
    sd.g1[is.na(sd.g1)]=0## if all bg, there's no info about fg sd
    ## add a shrinkage for sd.g1
    sd.prior=squeezeVar(sd.g1^2,n.g1-1)
    sd.g2=sqrt(sd.prior$var.post)
    #####  gene specific bg. Z_gi
    den.fg = den.bg = NA*Y
    for(i in 1:ncol(Y)){
        den.bg[,i]=dZinf.pois(Y[,i], par1[1,i], par1[2,i])
        den.fg[,i]=dLNP2(x=Y[,i], mu=mu.g1, sigma=sd.g2, l=L[i])
    }
    Z.fg=sweep(den.fg,2,1-pi0.hat,FUN="*")
    Z.bg=sweep(den.bg,2,pi0.hat,FUN="*")
    post.Z=Z.fg/(Z.fg+Z.bg)
    post.Z[is.na(post.Z)] <- 1

    ### if I shrink mu.g
    den.fg2 = NA*Y
    for (i in 1:ncol(Y)){
        den.fg2[,i]= dLNP2(x=Y[,i], mu=mu.g2, sigma=sd.g2, l=L[i])
    }
    Z.fg2=sweep(den.fg2,2,1-pi0.hat,FUN="*")
    post.Z2=Z.fg2/(Z.fg2+Z.bg)
    post.Z2[is.na(post.Z2)] <- 1

    pi.g = rowMeans(post.Z2)
    est_rslt = list(exprs = Y, pi.g = pi.g, p0 = par1[1,], lambda = par1[2,],
                    mu = mu.g2, sd = sd.g2, sf = L)
    return(est_rslt)
}

































