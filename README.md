<center> <h2> POWSC: A computational tool for power analysis in scRNA-seq </h2> </center>

power is a R package designed scRNA-seq with a wild range of usage. It can play three roles: parameter estimator, data simulator, power assessor. All the copyrights and explaned by Kenong Su <kenong.su@emroy.edu> and Dr. Wu's lab <http://www.haowulab.org>.

### 1. Software Installation
```
library(devtools)
install_github("suke18/POWSC", build_vignettes=TRUE)
```

### 2. Code Snippet
```
library(POWSC)
data("es.mef.eset")
eset = es.mef.eset[, pData(es.mef.eset)$celltype == "MEF"]
est_Paras = Eset2Phase(exprs(eset))
sim_size = c(100, 200, 500, 800, 1000)
################################################################
##### This is for the scenairo of pairwise cell comparison #####
################################################################
pow_rslt = runPOWSC(sim_size = sim_size, est_Paras = est_Paras, DE_Method = "MAST", Cell_Type = "PW")
plot(pow_rslt, Phase = "II")
summary(pow_rslt, Phase = "II")
```
