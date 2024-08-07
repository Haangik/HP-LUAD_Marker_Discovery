library(parallel)
library(pbapply)
library(dplyr)

##### Make a cluster for a parallel computing (optional)
##### Multiple core computing system with more than 100 cores is strongly recommended
##### It may take very long time if you run these codes with single core.
if(detectCores()>1){
  clu <- makePSOCKcluster(detectCores())  # You can modify this cluster object(clu) according to your computing environment
}else{
  clu <- NULL
}

clusterEvalQ(clu, {
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

## Univariate survival analysis by patient groups (Figure 2)
source("src/Univariate_survival_analysis.R")

### Univariate survival analysis with randomly sampled patients (Figure S1)
source("src/Patient_sampling_analysis.R")

## Risk prediction performance evaluation (Figure 3, Figure 4, and Table 2)
source("src/Performance_evaluation.R")

## Module enrichment analysis (Figure 5)
source("src/Module_enrichment_analysis.R")
