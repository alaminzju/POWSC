<center> <h2> POWSC: A computational tool for power analysis in scRNA-seq </h2> </center>

-------------------
`POWSC` is a R package designed for scRNA-seq with a wild range of usage. It can play three roles: parameter estimator, data simulator, and **power assessor**. As the parameter estimator, POWSC accurately captures the key characters to mimic the real expression data. As the data simulator, POWSC generates sythetic data based on a rigorous simulation mechanism incluidng zero expression values. As the power assessor, POWSC performs comprehensive power analysis and reports the stratified targeted power for two forms of DE genes. All the copyrights are explaned by Kenong Su <kenong.su@emroy.edu> and Dr. Wu's lab <http://www.haowulab.org>.
![Figure1](/assets/Figure1_3d40mf2bq.png)

### 1. Software Installation
```
library(devtools)
install_github("suke18/POWSC")
R CMD INSTALL POWSC_0.1.0.tar.gz
```

### 2. Code Snippet
**(1). parameter estimation for one cell type case**
```r
library(POWSC)
data("es_mef_sce")
sce = es_mef_sce[, colData(es_mef_sce)$cellTypes == "fibro"]
est_Paras = Est2Phase(sce) # Alternatively, you can load ("example_Paras"). The estimation is stored in est_Paras.
```

**(2). the first scenairo of two-group comparison**
```r
sim_size = c(50, 100, 200, 500, 1000)
pow_rslt = runPOWSC(sim_size = sim_size, est_Paras = est_Paras, DE_Method = "MAST", Cell_Type = "PW") # This might take a while (5 minutes).
plot(pow_rslt, Phase = "II") # Alternatively, using Phase = I
summary(pow_rslt, Phase = "II")
```
**(3). the second scenairo of multi-group comparisons**
```r
cell_per = c(0.1, 0.2, 0.3, 0.4)
powAll = vector(mode = "list", length = length(sim_size))
for (tmp_size in sim_size){
    tmpAll = runPOWSC(sim_size = tmp_size, est_Paras = estParas, DE_Method = "MAST", Cell_Type = "Multi", multi_Prob = cell_per)
    powAll[[toString(tmp_size)]] = tmpAll
}
```
**(4). A real data example**
```r
data("gbm_sce")
id = which(colData(gbm_sce)$patIds == "MGH26")
pat_sce = gbm_sce[, id]
sc3_rslt = sc3Mix(assays(pat_sce)$counts)
```
