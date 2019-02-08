<center> <h2> POWSC: A computational tool for power analysis in scRNA-seq </h2> </center>

POWSC is a R package designed scRNA-seq with a wild range of usage. It can play three roles: parameter estimator, data simulator, power assessor. All the copyrights and explaned by Kenong Su <kenong.su@emroy.edu> and Dr. Wu's lab <http://www.haowulab.org>.

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
data("es_mef_sce")
sce = es_mef_sce[, colData(es_mef_sce)$cellTypes == "MEF"]
est_Paras = Eset2Phase(sce)
sim_size = c(100, 200, 500, 800, 1000)
pow_rslt = runPOWSC(sim_size = sim_size, est_Paras = est_Paras, DE_Method = "MAST", Cell_Type = "PW")
plot(pow_rslt, Phase = "II")
summary(pow_rslt, Phase = "II")

################################################################
##### This is for the scenairo of multiple cell conditions #####
################################################################
cell_per = c(0.1, 0.2, 0.3, 0.4)
powAll = vector(mode = "list", length = length(sim_size))
for (tmp_size in sim_size){
    tmpAll = runPOWSC(sim_size = tmp_size, est_Paras = estParas, DE_Method = "MAST", Cell_Type = "Multi",
                      multi_Prob = cell_per)
    powAll[[toString(tmp_size)]] = tmpAll
}

################################################################
######### This is a real example from of scenairo two ##########
################################################################
data("gbm_sce")
id = which(colData(gbm_sce)$patIds == "MGH26")
pat_sce = gbm_sce[, id]
sc3_rslt = sc3Mix(assays(pat_sce)$counts)
est_Paras = Eset2Phase(pat_sce)
```
