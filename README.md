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
################################################################
##### This is for the scenairo of pairwise cell comparison #####
################################################################
data("es.mef.eset")
eset = es.mef.eset[, pData(es.mef.eset)$celltype == "MEF"]
est_Paras = Eset2Phase(exprs(eset))
sim_size = c(100, 200, 500, 800, 1000)
pow_rslt = runPOWSC(sim_size = sim_size, est_Paras = est_Paras, DE_Method = "MAST", Cell_Type = "PW")
plot(pow_rslt, Phase = "II")
summary(pow_rslt, Phase = "II")

################################################################
##### This is for the scenairo of multiple cell conditions #####
################################################################
cell_per = c(0.1, 0.2, 0.3, 0.4)
data("gbm.eset")
id = which(pData(gbm.eset)$patid == "patient id: MGH26")
mat = exprs(gbm.eset[, id])
estParas = Eset2Phase(mat)
powAll = vector(mode = "list", length = length(sim_size))
for (tmp_size in sim_size){
    tmpAll = runPOWSC(sim_size = tmp_size, est_Paras = estParas, DE_Method = "MAST", Cell_Type = "Multi",
                      multi_Prob = cell_per)
    powAll[[toString(tmp_size)]] = tmpAll
}
```
