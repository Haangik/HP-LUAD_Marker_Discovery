library(parallel)
library(pbapply)
library(dplyr)

##### Make a cluster for a parallel computing (optional)
##### Multiple core computing system with more than 100 cores is strongly recommended
##### It may take very long time if you run these codes with single core.
if(detectCores()>1){
  clu <- makePSOCKcluster(detectCores())  
}else{
  clu <- NULL
}

clusterEvalQ(clu, {
  setwd("/home/users/mogul012/data/Thesis")
  require(dplyr)
  require(effectsize)
  require(pROC)
  require(pwr)
  require(reshape2)
  require(expm)
  require(MetaIntegrator)
  require(meta)
  require(survival)
  require(metap)
  require(org.Hs.eg.db)
  NULL
})

## Data preparation
source("data/Patient_data_prep.R")
source("data/Functional_Module_prep.R")
source("data/Previous_published_markers.R")

## Statistical analysis
clusterExport(clu, "expression.data.list")
clusterExport(clu, "gene_dictionary")
clusterExport(clu, "clinical.annotation.list")
clusterExport(clu, c("na.terms", "dataset.name"))

egfr.datasets<-c('TCGA-LUAD', 'GSE11969', 'GSE13213', 'GSE26939', 'GSE31210', 'GSE72094')
egfr.dataset.ind<-match(egfr.datasets, dataset.name)

## Module enrichment analysis

## Performance evaluation
