######## POWSC ######## 
runPOWSC = function(sim_size = c(50, 100, 200, 800, 1000), per_DE = 0.05, est_Paras, DE_Method = c("MAST", "SC2P"), 
                    Cell_Type = c("PW", "Multi"), multi_Prob = NULL, 
                    alpha = 0.1, disc_delta = 0.1, cont_delta = 0.5){
    DE_Method = match.arg(DE_Method)
    Cell_Type = match.arg(Cell_Type)
    pow_rslt = NULL
    # Pair-wise comparison 
    if (Cell_Type == "PW"){
        for (tmp_size in sim_size){
            sim_Data = Simulate2SCE(n = tmp_size, estParas1 = est_Paras, estParas2 = est_Paras)
            DE_rslt = runDE(sim_Data$sce, DE_Method = DE_Method)
            pow1_rslt = Power_Disc(DErslt = DE_rslt, simData = sim_Data, alpha = alpha, delta = disc_delta)
            pow2_rslt = Power_Cont(DErslt = DE_rslt, simData = sim_Data, alpha = alpha, delta = cont_delta)
            tmp_rslt = list(pow1 = pow1_rslt, pow2 = pow2_rslt)
            pow_rslt[[toString(tmp_size)]] = tmp_rslt
        }
    }else{ # Multiple Cell Types
        if (missing(multi_Prob)){
            stop("Please Estimate the Cell Proportions by SC3 or Seurat")
        }
        for (tmp_size in sim_size){
            pow1_rslt = pow2_rslt = NULL
            sim_Data = SimulateMultiSCEs(n = tmp_size, estParas = est_Paras, multiProb = multi_Prob)
            for (comp in names(sim_Data)){
                tmp_DE = runDE(sim_Data[[comp]]$sce, DE_Method = DE_Method)
                pow1_rslt[[comp]] = Power_Disc(DErslt = tmp_DE, simData = sim_Data[[comp]], alpha = alpha, delta = disc_delta)
                pow2_rslt[[comp]] = Power_Cont(DErslt = tmp_DE, simData = sim_Data[[comp]], alpha = alpha, delta = cont_delta)
            }
            tmp_rslt = list(pow1 = pow1_rslt, pow2 = pow2_rslt)
            pow_rslt[[toString(tmp_size)]] = tmp_rslt
        }
    }
    class(pow_rslt) = "POWSC"
    return(pow_rslt)
}

######## Plot ########
plot.POWSC = function(POWSCobj, Phase = c("I", "II")){
    pow1_mat = do.call(rbind, lapply(POWSCobj, function(x) x$pow1$power))
    pow2_mat = do.call(rbind, lapply(POWSCobj, function(x) x$pow2$power))
    nrep = rownames(pow1_mat)
    nm1 = names(POWSCobj[[1]][[1]][[1]])
    nm2 = names(POWSCobj[[1]][[2]][[1]])
    pow1 = data.frame(Strata = rep(nm1, each = nrow(pow1_mat)), Power = as.vector(pow1_mat), Reps = rep(nrep, ncol(pow1_mat)))
    pow2 = data.frame(Strata = rep(nm2, each = nrow(pow2_mat)), Power = as.vector(pow2_mat), Reps = rep(nrep, ncol(pow2_mat)))
    pow1$Strata = factor(pow1$Strata, levels = nm1); pow2$Strata = factor(pow2$Strata, levels = nm2)
    Ph = match.arg(Phase)
    if (Ph == "I"){
        pow = pow1
        tit = ggtitle("DE from Phase I")
        tmpxlab = xlab("strata of zero ratios across samples")
    }else{
        pow = pow2
        tit = ggtitle("DE from Phase II")
        tmpxlab = xlab("strata of average reads")
    }
    gp = ggplot(pow, aes(x=Strata, y=Power, group=Reps, color=Reps)) +
        geom_line(aes(color=Reps, linetype = Reps), size = 0.6)+
        geom_point(aes(color=Reps, shape = Reps), size = 1.5) +
        scale_colour_brewer(palette = "Set1") +
        tmpxlab + ylab("power") + theme_classic() +
        scale_y_continuous(breaks=seq(0,1,by=0.1),limits = c(0, 1), labels=seq(0,1,by=0.1))+
        theme(legend.position = c(0.5, 0.1), legend.direction = "horizontal") +
        theme(plot.title = element_text(size=16, face="bold"))+ tit
    plot(gp)
}

######## Summary Table ########
summary.POWSC = function(POWSCobj, Phase = c("I", "II")){
    pow1_mat = round(do.call(rbind, lapply(POWSCobj, function(x) x$pow1$power)), 4)
    pow2_mat = round(do.call(rbind, lapply(POWSCobj, function(x) x$pow2$power)), 4)
    nm1 = names(POWSCobj[[1]][[1]][[1]])
    nm2 = names(POWSCobj[[1]][[2]][[1]])
    colnames(pow1_mat) = nm1
    colnames(pow2_mat) = nm2
    Ph = match.arg(Phase)
    if (Ph == "I"){
        tab = pow1_mat
    }else{
        tab = pow2_mat
    }
    return(tab)
}













