<center> <h2> POWSC: A computational tool for power analysis in scRNA-seq </h2> </center>

-------------------
**_POWSC_** is a R package designed for scRNA-seq with a wild range of usage. It can play three roles: **parameter estimator**, **data simulator**, and **power assessor**. As the parameter estimator, POWSC accurately captures the key characters (`Fig.B`) for any specific cell type from a given real expression data (`Fig.A`). As the data simulator, POWSC generates sythetic data (`Fig.C`) based on a rigorous simulation mechanism incluidng zero expression values. As the power assessor, POWSC performs comprehensive power analysis and reports the stratified targeted power (`Fig.D`) for two forms of DE genes. A schemetic overview of the aglorithm is shown in (`Fig.E`). All the copyrights are explaned by Kenong Su <kenong.su@emroy.edu> and Dr. Wu's lab <http://www.haowulab.org>.
![workflow](/assets/workflow.png)

### 1. Software Installation
```
library(devtools)
install_github("suke18/POWSC")
R CMD INSTALL POWSC_0.1.0.tar.gz
```

### 2. Code Snippets
**(1). parameter estimation for one cell type case**
```r
library(POWSC)
data("es_mef_sce")
sce = es_mef_sce[, colData(es_mef_sce)$cellTypes == "fibro"]
set.seed(1)
# For demonstration purpose, we extract a small dataset here:
sce_small = sce[sample(1:nrow(sce), 5000),]
est_Paras = Est2Phase(sce_small)
```

**(2). the first scenairo of two-group comparison**
```r
sim_size = c(100, 400, 1000) # A numeric vector
pow_rslt = runPOWSC(sim_size = sim_size, est_Paras = est_Paras,per_DE=0.05, DE_Method = "MAST", Cell_Type = "PW")
plot(pow_rslt, Form="II", Cell_Type = "PW") # Alternatively, we can use Form="I"
summary(pow_rslt, Form="II", Cell_Type = "PW")
```
**(3). the second scenairo of multi-group comparisons**
```r
sim_size = 1000
cell_per = c(0.2, 0.3, 0.5)
pow_rslt = runPOWSC(sim_size = sim_size, est_Paras = est_Paras,per_DE=0.05, DE_Method = "MAST",multi_Prob = cell_per, Cell_Type = "Multi")
plot(pow_rslt, Form="II", Cell_Type = "Multi")
summary(pow_rslt, Form="II", Cell_Type = "Multi")
```
