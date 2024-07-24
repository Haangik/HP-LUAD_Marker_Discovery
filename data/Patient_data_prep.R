library(org.Hs.eg.db)
library(SummarizedExperiment)
library(tidyverse)
library(TCGAbiolinks)
library(GEOquery)

# Terms representing elements not available were defined in advance of the analysis
na.terms<-c(NA, "[Discrepancy]", "[Not Applicable]", "[Not Available]", "[Not Evaluated]", "[Unknown]", "MX", "TX", "NX")

alive.days.thres=0*365
alive.month.thres=alive.days.thres/30.41667

death.days.thres=100000*365
death.month.thres=death.days.thres/30.41667


# 1. TCGA data
target_cancer <- "TCGA-LUAD"

###### Clinical data
query<-GDCquery(
  project = target_cancer,
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  file.type = ".txt",
)

GDCdownload(
  query = query
)

clinical<-GDCprepare(query)
clinical<-lapply(clinical, FUN=function(x){
  return(x[-c(1,2),])
})

####### Mutation data
query <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation", 
)
GDCdownload(query)
maf <- GDCprepare(query)


######## Expression data
query <- GDCquery(
  project = target_cancer,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)

clinical_df<-data.frame(data@colData@listData[c('barcode', 'sample_type', 'ajcc_pathologic_stage', 'age_at_index',
                                                "days_to_last_follow_up", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m", 
                                                "ajcc_staging_system_edition", 'race', 'gender', "vital_status", "days_to_death", "paper_expression_subtype")])

colnames(clinical_df)[4]<-"age"

clinical_df$survival_time<-sapply(c(1:length(clinical_df$vital_status)), FUN=function(x){
  t=setdiff(c(clinical_df$days_to_last_follow_up[x], clinical_df$days_to_death[x]), na.terms)
  if(length(t)==0){
    return("[Not Available]")
  }
  return(max(as.numeric(t)))
})

clinical_df$relapse<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  if(is.na(match(barcode.cut, clinical[[2]]$bcr_patient_barcode))){
    return("NO")
  }else{
    return("YES")
  }
})

clinical_df$race<-sapply(clinical_df$race, FUN=function(x){
  if(x=="not reported")
    return("[Not Available]")
  else
    return(x)
})

clinical_df$adjuvant_chemo<-sapply(data@colData@listData$treatments, FUN=function(x){
  res<-x %>% filter(treatment_type=="Pharmaceutical Therapy, NOS") %>% dplyr::select(treatment_or_therapy) %>% pull
  if(res=="not reported"){
    res<-"[Not Available]"
  }
  return(res)
})

clinical_df$adjuvant_radiation<-sapply(data@colData@listData$treatments, FUN=function(x){
  res<-x %>% filter(treatment_type=="Radiation Therapy, NOS") %>% dplyr::select(treatment_or_therapy) %>% pull
  if(res=="not reported"){
    res<-"[Not Available]"
  }
  return(res)
})

tumor.ind<-which(clinical_df$sample_type=="Primary Tumor")

clinical_df<-clinical_df[tumor.ind,]

clinical_df$barcode.short<-sapply(clinical_df$barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

# smoking status
clinical_df$smoking_status<-sapply(clinical_df$barcode.short, FUN=function(barcode){
  ind<-match(barcode, clinical$clinical_patient_luad$bcr_patient_barcode)
  smoking.indicator<-clinical$clinical_patient_luad$tobacco_smoking_history_indicator[ind]
  if(smoking.indicator == "1"){
    res<-"Never smoker"
  }else if (smoking.indicator %in% c("2", "3", "4", "5")){
    res<-"Ever smoker"
  }else{
    res<-"[Not Available]"
  }
  return(res)
})

# Predominant subtype
clinical_df$predominant_subtype<-sapply(clinical_df$barcode.short, FUN=function(barcode){
  ind<-match(barcode, clinical$clinical_patient_luad$bcr_patient_barcode)
  x<-clinical$clinical_patient_luad$histologic_diagnosis...68[ind]
  if(x=="Lung Acinar Adenocarcinoma"){
    "Acinar"
  }else if(x=="Lung Micropapillary Adenocarcinoma"){
    "Micropapillary"
  }else if(x %in% c("Lung Bronchioloalveolar Carcinoma Mucinous", "Lung Mucinous Adenocarcinoma", "Mucinous (Colloid) Carcinoma")){
    "Mucinous"
  }else if(x=="Lung Papillary Adenocarcinoma"){
    "Papillary"
  }else{
    "[Not Available]"
  }
  
})

# EGFR mutation
mut.egfr<-maf %>% filter(Hugo_Symbol == "EGFR")
mut.egfr$barcode.short<-sapply(mut.egfr$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df$egfr_mut<-sapply(clinical_df$barcode.short, FUN=function(x){
  ind<-match(x, mut.egfr$barcode.short)
  if(is.na(ind))
    res<-"EGFR-WT"
  else
    res<-"EGFR-Mut"
})

# KRAS mutation
mut.kras<-maf %>% filter(Hugo_Symbol == "KRAS")
mut.kras$barcode.short<-sapply(mut.kras$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df$kras_mut<-sapply(clinical_df$barcode.short, FUN=function(x){
  ind<-match(x, mut.kras$barcode.short)
  if(is.na(ind))
    res<-"KRAS-WT"
  else
    res<-"KRAS-Mut"
})

# p53 mutation
mut.tp53<-maf %>% filter(Hugo_Symbol == "TP53")
mut.tp53$barcode.short<-sapply(mut.tp53$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df$p53_mut<-sapply(clinical_df$barcode.short, FUN=function(x){
  ind<-match(x, mut.tp53$barcode.short)
  if(is.na(ind))
    res<-"p53-WT"
  else
    res<-"p53-Mut"
})

# STK11 mutation
mut.stk11<-maf %>% filter(Hugo_Symbol == "STK11")
mut.stk11$barcode.short<-sapply(mut.stk11$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df$stk11_mut<-sapply(clinical_df$barcode.short, FUN=function(x){
  ind<-match(x, mut.stk11$barcode.short)
  if(is.na(ind))
    res<-"STK11-WT"
  else
    res<-"STK11-Mut"
})


data.luad <- assay(data, "fpkm_uq_unstrand")
data.luad.row.info<-rowRanges(data)
protein_coding_gene.ind<-which(data.luad.row.info@elementMetadata@listData$gene_type=="protein_coding")

data.luad.protein_coding<-data.luad[protein_coding_gene.ind, tumor.ind]

expression_df<-t(data.luad.protein_coding)
colnames(expression_df)<-row.names(data.luad.protein_coding)
expression_df<-log(expression_df+1, base=2)

gene_info<-data.frame(gene_id=data.luad.row.info@elementMetadata@listData$gene_id, 
                      gene_symbol=data.luad.row.info@elementMetadata@listData$gene_name) %>% filter(gene_id %in% colnames(expression_df))

gene_info$ENTREZID<-pbsapply(gene_info$gene_symbol, FUN=function(x){
  tryCatch(AnnotationDbi::select(org.Hs.eg.db, 
                                 keys=x, 
                                 columns=c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")$ENTREZID[1], 
           error= function(e){ 
             return(NA)
           })
}, cl=clu)

clusterExport(clu, "gene_info")


gene_info$ENTREZID<-pbsapply(c(1:length(gene_info$ENTREZID)), FUN=function(k){
  if(!is.na(gene_info$ENTREZID[k])){
    return(gene_info$ENTREZID[k])
  }else{
    tryCatch(AnnotationDbi::select(org.Hs.eg.db, 
                                   keys=unlist(strsplit(gene_info$gene_id[k], ".", fixed=T))[1], 
                                   columns=c("ENSEMBL", "ENTREZID"), keytype = "ENSEMBL")$ENTREZID[1], 
             error= function(e){ 
               return(NA)
             })
  }
}, cl=clu)


clinical_TCGA<-clinical_df
expression_TCGA<-expression_df

rm(clinical_df)
rm(expression_df)
rm(data.luad); rm(data.luad.protein_coding); rm(data.Platform); rm(data.luad.row.info)
rm(maf);rm(mut.egfr);rm(mut.kras);rm(mut.stk11);rm(mut.tp53);rm(query)


# 3. GSE11969
data <- getGEO(filename="GEO/GSE11969_series_matrix.txt.gz")
clindata <- data@phenoData@data

clinical_GSE11969<-data.frame(Accession=clindata$geo_accession,
                              age=clindata$`Age:ch1`,
                              histology=sapply(clindata$characteristics_ch1.4, FUN=function(x){
                                if(x=="Histology: AD"){
                                  "ADC"
                                }else if(x=="Histology: SQ"){
                                  "SCC"
                                }else if(x=="Histology: Normal lung tissue"){
                                  "Normal"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              gender=sapply(clindata$characteristics_ch1.3, FUN=function(x){
                                if(x=="Sex: F"){
                                  "female"
                                }else if(x=="Sex: M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_t=sapply(str_sub(clindata$`pTNM:ch1`,1,2), FUN=function(x){
                                if(x %in% c("T1","T2", "T3", "T4", "TX")){
                                  x
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_n=sapply(str_sub(clindata$`pTNM:ch1`,3,4), FUN=function(x){
                                if(x %in% c("NX", "N0", "N1", "N2", "N3")){
                                  x
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              # thie dataset has ajcc_pathologic_m label but only contains M0
                              ajcc_pathologic_stage=sapply(clindata$`pStage:ch1`, FUN=function(x){
                                if(x=="Stage: 1"){
                                  "Stage I"
                                }else if(x=="IA"){
                                  "Stage IA"
                                }else if(x=="IB"){
                                  "Stage IB"
                                }else if(x=="IIA"){
                                  "Stage IIA"
                                }else if(x=="IIB"){
                                  "Stage IIB"
                                }else if(x=="IIIA"){
                                  "Stage IIIA"
                                }else if(x=="IIIB"){
                                  "Stage IIIB"
                                }else if(x=="IV"){
                                  "Stage IV"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              egfr_mut=sapply(clindata$`EGFR status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "EGFR-WT"
                                }else if(x=="Mut"){
                                  "EGFR-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              kras_mut=sapply(clindata$`K-ras Status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "KRAS-WT"
                                }else if (x=="Mut"){
                                  "KRAS-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              p53_mut=sapply(clindata$`p53 Status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "p53-WT"
                                }else if(x=="Mut"){
                                  "p53-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              grade=sapply(clindata$`Grade:ch1`, FUN=function(x){
                                if(x=="Moderately"){
                                  "Moderately"
                                }else if(x=="Poorly"){
                                  "Poorly"
                                }else if (x=="Well"){
                                  "Well"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$`Survival (days):ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              vital_status=sapply(clindata$`Status:ch1`, FUN=function(x){
                                if(x=="Alive"){
                                  "Alive"
                                }else if(x=="Dead"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              relapse=sapply(clindata$`Evidence of relapse:ch1`, FUN=function(x){
                                if(is.na(x))
                                  return("[Not Available]")
                                if(x=="N"){
                                  "NO"
                                }else if(x=="Y"){
                                  "YES"
                                }else
                                  "[Not Available]"
                              })
)

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data


# Mapping Probe ID to Symbol
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Symbol=unlist(strsplit(data.Platform$`Gene symbol`[k], " ",fixed = T))
  if(length(Symbol)!=0){
    df<-data.frame(ID, Symbol)  
  }else{
    df<-data.frame(ID=ID, Symbol=NA)
  }
  return(df)
})) %>% filter(!is.na(Symbol))

data.MA<-data.MA[mapping.table$ID,]

symbol.to.id<-pblapply(unique(mapping.table$Symbol), FUN=function(geneid){
  indices=which(mapping.table$Symbol %in% geneid)
})
names(symbol.to.id)<-unique(mapping.table$Symbol)

data.MA<-do.call(rbind, pblapply(c(1:length(symbol.to.id)), FUN=function(k){
  if(length(symbol.to.id[[k]])>1){
    return(colMeans(data.MA[symbol.to.id[[k]],]))  
  }else{
    return(data.MA[symbol.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(symbol.to.id)


expression_GSE11969<-t(data.MA)
adc.ind<-which(clinical_GSE11969$histology=="ADC")
clinical_GSE11969<-clinical_GSE11969[adc.ind,]
expression_GSE11969<-expression_GSE11969[adc.ind, ]


# 4. GSE13213
data <- getGEO(filename="GEO/GSE13213_series_matrix.txt.gz")
clindata <- data@phenoData@data

clinical_GSE13213<-data.frame(Accession=clindata$geo_accession,
                              age=clindata$`Age:ch1`,
                              histology=sapply(clindata$characteristics_ch1.4, FUN=function(x){
                                if(x=="Histology: AD"){
                                  "ADC"
                                }else if(x=="Histology: SQ"){
                                  "SCC"
                                }else if(x=="Histology: Normal lung tissue"){
                                  "Normal"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              gender=sapply(clindata$characteristics_ch1.3, FUN=function(x){
                                if(x=="Sex: F"){
                                  "female"
                                }else if(x=="Sex: M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_t=sapply(str_sub(clindata$`TNM (Pathological):ch1`,1,2), FUN=function(x){
                                if(x %in% c("T1","T2", "T3", "T4", "TX")){
                                  x
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_n=sapply(str_sub(clindata$`TNM (Pathological):ch1`,3,4), FUN=function(x){
                                if(x %in% c("NX", "N0", "N1", "N2", "N3")){
                                  x
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(clindata$`Stage (Pathological ):ch1`, FUN=function(x){
                                if(x=="IA"){
                                  "Stage IA"
                                }else if(x=="IB"){
                                  "Stage IB"
                                }else if(x=="IIA"){
                                  "Stage IIA"
                                }else if(x=="IIB"){
                                  "Stage IIB"
                                }else if(x=="IIIA"){
                                  "Stage IIIA"
                                }else if(x=="IIIB"){
                                  "Stage IIIB"
                                }else if(x=="IV"){
                                  "Stage IV"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              egfr_mut=sapply(clindata$`EGFR status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "EGFR-WT"
                                }else if(x=="Mut"){
                                  "EGFR-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              kras_mut=sapply(clindata$`K-ras Status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "KRAS-WT"
                                }else if (x=="Mut"){
                                  "KRAS-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              p53_mut=sapply(clindata$`p53 Status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "p53-WT"
                                }else if(x=="Mut"){
                                  "p53-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$`Survival (days):ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              vital_status=sapply(clindata$`Status:ch1`, FUN=function(x){
                                if(x=="Alive"){
                                  "Alive"
                                }else if(x=="Dead"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              relapse=sapply(clindata$`Evidence of relapse:ch1`, FUN=function(x){
                                if(is.na(x))
                                  return("[Not Available]")
                                if(x=="N"){
                                  "NO"
                                }else if(x=="Y"){
                                  "YES"
                                }else
                                  "[Not Available]"
                              })
)

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Removing probes which contains NA values
na.ind<-which(apply(data.MA,MARGIN=1, FUN=function(x){
  length(which(is.na(x)))
})>0)
data.MA<-data.MA[-na.ind,]
data.Platform<-data.Platform[-na.ind,]

# Mapping Probe ID to Symbol
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$GENE[k]
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]


entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(na.omit(data.MA[entrez.to.id[[k]],])))
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE13213<-data.MA
adc.ind<-which(clinical_GSE13213$histology=="ADC")
clinical_GSE13213<-clinical_GSE13213[adc.ind,]
expression_GSE13213<-t(expression_GSE13213[, adc.ind])



# 7. GSE26939
data <- getGEO(filename="GEO/GSE26939_series_matrix.txt.gz")
clindata <- data@phenoData@data

clinical_GSE26939<-data.frame(Accession=clindata$geo_accession,
                              age=clindata$`age (90=greater than or equal to 90):ch1`,
                              gender=sapply(c(1:dim(clindata)[1]), FUN=function(k){
                                x=setdiff(c(clindata$`sex:ch1`[k], clindata$`Sex:ch1`[k]), NA)
                                if(x=="F"){
                                  "female"
                                }else if(x=="M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              smoking_status=sapply(clindata$`smoking_status(0=nonsmoker,1=smoker):ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="Smoker"){
                                  "Ever smoker"
                                }else if(x=="NeverSmoker"){
                                  "Never smoker"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(c(1:dim(clindata)[1]), FUN=function(k){
                                x=setdiff(c(clindata$`stage:ch1`[k], clindata$`Stage:ch1`[k]), NA)
                                if(length(x)==0){
                                  return("[Not Available]")
                                }
                                if(x=="IA"){
                                  "Stage IA"
                                }else if(x=="IB"){
                                  "Stage IB"
                                }else if(x=="IIA"){
                                  "Stage IIA"
                                }else if(x=="IIB"){
                                  "Stage IIB"
                                }else if(x=="IIIA"){
                                  "Stage IIIA"
                                }else if(x=="IIIB"){
                                  "Stage IIIB"
                                }else if(x=="IV"){
                                  "Stage IV"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              egfr_mut=sapply(clindata$`egfr (0='wt',1='mutated'):ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="1"){
                                  "EGFR-Mut"
                                }else if(x=="0"){
                                  "EGFR-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              kras_mut=sapply(clindata$`kras (0='wt',1='mutated'):ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="1"){
                                  "KRAS-Mut"
                                }else if(x=="0"){
                                  "KRAS-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              stk11_mut=sapply(clindata$`stk11 (0='wt',1='mutated'):ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="1"){
                                  "STK11-Mut"
                                }else if(x=="0"){
                                  "STK11-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              predimonant_subtype=sapply(clindata$`subtype:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }else{
                                  return(x)
                                }
                              }),
                              survival_time=sapply(clindata$`survival_months:ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              vital_status=sapply(clindata$`survival_status(0='alive',1='dead'):ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="0"){
                                  "Alive"
                                }else if(x=="1"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              })
)



data.MA=data@assayData$exprs
data.Platform<-data@featureData@data


# Mapping Probe ID to Symbol
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Symbol=unlist(strsplit(data.Platform$ORF[k], " ",fixed = T))
  if(length(Symbol)!=0){
    df<-data.frame(ID, Symbol)  
  }else{
    df<-data.frame(ID=ID, Symbol=NA)
  }
  return(df)
})) %>% filter(!is.na(Symbol))

data.MA<-data.MA[match(mapping.table$ID, rownames(data.MA)),]

symbol.to.id<-pblapply(unique(mapping.table$Symbol), FUN=function(geneid){
  indices=which(mapping.table$Symbol %in% geneid)
})
names(symbol.to.id)<-unique(mapping.table$Symbol)

data.MA<-do.call(rbind, pblapply(c(1:length(symbol.to.id)), FUN=function(k){
  if(length(symbol.to.id[[k]])>1){
    return(colMeans(data.MA[symbol.to.id[[k]],]))  
  }else{
    return(data.MA[symbol.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(symbol.to.id)

expression_GSE26939<-t(data.MA)


# 10. GSE31210
data <- getGEO(filename="GEO/GSE31210_series_matrix.txt.gz")
clindata <- data@phenoData@data

clinical_GSE31210<-data.frame(Accession=clindata$geo_accession,
                              age=clindata$`age (years):ch1`,
                              gender=sapply(clindata$`gender:ch1`, FUN=function(x){
                                if(x=="female"){
                                  "female"
                                }else if(x=="male"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              histology=sapply(clindata$`tissue:ch1`, FUN=function(x){
                                if(x=="primary lung tumor"){
                                  "ADC"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(clindata$`pathological stage:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x == "IA"){
                                  "Stage IA"
                                }else if(x == "IB"){
                                  "Stage IB"
                                }else if(x == "II"){
                                  "Stage II"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              egfr_mut=sapply(clindata$`gene alteration status:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="EGFR mutation +"){
                                  "EGFR-Mut"
                                }else{
                                  "EGFR-WT"
                                }
                              }),
                              kras_mut=sapply(clindata$`gene alteration status:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="KRAS mutation +"){
                                  "KRAS-Mut"
                                }else{
                                  "KRAS-WT"
                                }
                              }),
                              survival_time=sapply(clindata$`days before death/censor:ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              vital_status=sapply(clindata$`death:ch1` , FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="alive"){
                                  "Alive"
                                }else if(x=="dead"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              disease_free_survival=sapply(clindata$`days before relapse/censor:ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              relapse=sapply(clindata$`relapse:ch1`, FUN=function(x){
                                if(is.na(x))
                                  return("[Not Available]")
                                if(x=="not relapsed"){
                                  "NO"
                                }else if(x=="relapsed"){
                                  "YES"
                                }else
                                  "[Not Available]"
                              })
)

data.MA=data@assayData$exprs %>% log(2)
data.Platform<-data@featureData@data
# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
entrez.to.id<-lapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE31210<-t(data.MA)
adc.ind<-which(clinical_GSE31210$histology=="ADC")
clinical_GSE31210<-clinical_GSE31210[adc.ind,]
expression_GSE31210<-expression_GSE31210[adc.ind,]


# 17. GSE72094
data <- getGEO(filename="GEO/GSE72094_series_matrix.txt.gz")
clindata <- data@phenoData@data
clinical_GSE72094<-data.frame(Accession=clindata$geo_accession,
                              age=clindata$`age_at_diagnosis:ch1`,
                              gender=sapply(clindata$`gender:ch1`, FUN=function(x){
                                if(x=="F"){
                                  "female"
                                }else if(x=="M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              race=sapply(clindata$`race:ch1`, FUN=function(x){
                                if(x %in% c("ASIAN INDIAN OR PAKISTANI", "VIETNAMESE", "THAI")){
                                  "asian"
                                }else if(x=="BLACK"){
                                  "black or african american"
                                }else if(x=="WHITE"){
                                  "white"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              smoking_status=sapply(clindata$`smoking_status:ch1`, FUN=function(x){
                                if(x== "Ever"){
                                  "Ever smoker"
                                }else if(x=="Never"){
                                  "Never smoker"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              
                              ajcc_pathologic_stage=sapply(clindata$characteristics_ch1.11, FUN=function(x){
                                if(x=="Stage: 1"){
                                  "Stage I"
                                }else if(x=="Stage: 1A"){
                                  "Stage IA"
                                }else if(x=="Stage: 1B"){
                                  "Stage IB"
                                }else if(x=="Stage: 2A"){
                                  "Stage IIA"
                                }else if(x=="Stage: 2B"){
                                  "Stage IIB"
                                }else if(x=="Stage: 3"){
                                  "Stage III"
                                }else if(x=="Stage: 3A"){
                                  "Stage IIIA"
                                }else if(x=="Stage: 3B"){
                                  "Stage IIIB"
                                }else if(x=="Stage: 4"){
                                  "Stage IV"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              egfr_mut=sapply(clindata$`egfr_status:ch1`, FUN=function(x){
                                if(x=="Mut"){
                                  "EGFR-Mut"
                                }else if(x=="WT"){
                                  "EGFR-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              kras_mut=sapply(clindata$`kras_status:ch1`, FUN=function(x){
                                if(x=="Mut"){
                                  "KRAS-Mut"
                                }else if(x=="WT"){
                                  "KRAS-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              p53_mut=sapply(clindata$`tp53_status:ch1`, FUN=function(x){
                                if(x=="Mut"){
                                  "p53-Mut"
                                }else if(x=="WT"){
                                  "p53-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              stk11_mut=sapply(clindata$`stk11_status:ch1`, FUN=function(x){
                                if(x=="Mut"){
                                  "STK11-Mut"
                                }else if(x=="WT"){
                                  "STK11-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$characteristics_ch1.10, FUN=function(x){
                                res<-as.numeric(gsub("survival_time_in_days: ","", x))
                                if(is.na(res))
                                  res<-"[Not Available]"
                                return(res)
                              }),
                              vital_status=sapply(clindata$characteristics_ch1.9, FUN=function(x){
                                res<-gsub("vital_status: ","", x)
                                if(res=="NA")
                                  res<-"[Not Available]"
                                return(res)
                              }))


data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$EntrezGeneID[k], " ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]

entrez.to.id<-lapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE72094<-t(data.MA)

rm(data);rm(data.MA); rm(entrez.to.id);
rm(mapping.table); rm(adc.ind);
rm(alive.days.thres);rm(alive.month.thres);rm(death.days.thres);rm(death.month.thres)

################# gene dictionary generation
## checking validity of gene symbol of TCGA data
write.csv(gene_info$gene_symbol,"gene_symbols_tcga.csv", row.names=F, quote=F)
gene_symbols.tcga<-read.csv("hgnc-symbol-check-tcga.csv")

clusterEvalQ(clu, {
  gene_symbols.tcga<-read.csv("hgnc-symbol-check-tcga.csv")
  NULL
})
gene_info$gene_symbol<-pbsapply(c(1:length(gene_info$gene_symbol)), FUN=function(i){
  symbol.info<-gene_symbols.tcga %>% filter(Input %in% gene_info$gene_symbol[i])
  if(!is.na(match("Approved symbol", symbol.info$Match.type))){
    res<-symbol.info$Approved.symbol[match("Approved symbol", symbol.info$Match.type)]
  }else{
    if(dim(symbol.info)[1]==1){
      res<-symbol.info$Approved.symbol
    }else{
      res<-symbol.info$Approved.symbol[match("Previous symbol", symbol.info$Match.type)]
    }
  }
  return(res)
}, cl=clu)

clusterExport(clu, "gene_info")

gene_info$ENTREZID<-pbsapply(c(1:length(gene_info$gene_symbol)), FUN=function(i){
  if(!is.na(gene_info$ENTREZID[i])){
    return(gene_info$ENTREZID[i])
  }else{
    tryCatch(AnnotationDbi::select(org.Hs.eg.db, 
                                   keys=gene_info$gene_symbol[i], 
                                   columns=c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")$ENTREZID, 
             error= function(e){ 
               return(NA) 
             })
  }
}, cl=clu)


whole.gene.symbol<-Reduce(union, list(colnames(expression_GSE11969),colnames(expression_GSE26939)))
whole.gene.entrez<-Reduce(union, list(colnames(expression_GSE13213), colnames(expression_GSE31210), 
                                      colnames(expression_GSE72094)))

# gene ids with symbols
gene_ids<-AnnotationDbi::select(org.Hs.eg.db, 
                                keys=setdiff(whole.gene.entrez, gene_info$ENTREZID), 
                                columns=c("SYMBOL", "ENTREZID"), keytype = "ENTREZID")

for(i in match(intersect(gene_ids$SYMBOL, gene_info$gene_symbol), gene_info$gene_symbol)){
  gene_info$ENTREZID[i]<-gene_ids$ENTREZID[match(intersect(gene_ids$SYMBOL, gene_info$gene_symbol), gene_ids$SYMBOL)]
}

gene_ids<-gene_ids %>% filter(ENTREZID %in% setdiff(whole.gene.entrez, gene_info$ENTREZID))
gene_ids<-data.frame(gene_id=NA,
                     gene_symbol=gene_ids$SYMBOL,
                     ENTREZID=gene_ids$ENTREZID)

write.csv(na.omit(gene_ids$gene_symbol),"gene_symbols_entrezid.csv", row.names=F, quote=F)
gene_symbols.entrez<-read.csv("hgnc-symbol-check-entrez.csv")

clusterEvalQ(clu, {
  gene_symbols.entrez<-read.csv("hgnc-symbol-check-entrez.csv")
  NULL
})
clusterExport(clu, c("gene_ids", "gene_symbols.entrez"))

gene_ids$gene_symbol<-pbsapply(c(1:length(gene_ids$gene_symbol)), FUN=function(i){
  if(is.na(gene_ids$gene_symbol[i])){
    return(NA)
  }
  
  symbol.info<-gene_symbols.entrez %>% filter(Input %in% gene_ids$gene_symbol[i])
  if(!is.na(match("Approved symbol", symbol.info$Match.type))){
    res<-symbol.info$Approved.symbol[match("Approved symbol", symbol.info$Match.type)]
  }else{
    if(dim(symbol.info)[1]==1){
      res<-symbol.info$Approved.symbol
    }else{
      symbol.info<-symbol.info %>% filter(!Location %in% "mitochondria")
      res<-symbol.info$Approved.symbol[1]
    }
  }
  return(res)
}, cl=clu)


# gene symbol with entrez ids
gene_symbols<-AnnotationDbi::select(org.Hs.eg.db,
                                    keys=setdiff(setdiff(whole.gene.symbol, gene_info$gene_symbol), gene_ids$SYMBOL),
                                    columns=c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")

write.csv(setdiff(setdiff(whole.gene.symbol, gene_info$gene_symbol), gene_ids$SYMBOL),
          "gene_symbols.csv", row.names=F, quote=F)

### HUGO multi-symbol checker
symbol.dict<-read.csv("hgnc-symbol-check.csv")

for(i in 1:dim(gene_symbols)[1]){
  symbol.info<-symbol.dict %>% filter(Input %in% gene_symbols$SYMBOL[i])
  if(dim(symbol.info)[1]==0){
    next
  }
  
  if(symbol.info$Match.type[1]=="Unmatched"){
    next
  }
  
  if(dim(symbol.info)[1]==1){
    symbol.update<-symbol.info$Approved.symbol
  }else{
    if(!is.na(match("Approved symbol", symbol.info$Match.type))){
      symbol.update<-symbol.info$Approved.symbol[which(symbol.info$Match.type=="Approved symbol")]
    }else{
      symbol.update<-symbol.info$Approved.symbol[1]
    }
  }

  if(!is.na(match(gene_symbols$SYMBOL[i], colnames(expression_GSE26939)))){
    colnames(expression_GSE26939)[match(gene_symbols$SYMBOL[i], colnames(expression_GSE26939))]<-symbol.update
  }
}

whole.gene.symbol<-Reduce(union, list(colnames(expression_GSE11969), colnames(expression_GSE19188), colnames(expression_GSE26939)))
gene_symbols<-AnnotationDbi::select(org.Hs.eg.db,
                                    keys=setdiff(whole.gene.symbol, c(gene_info$gene_symbol, gene_ids$SYMBOL)),
                                    columns=c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")
gene_symbols<-gene_symbols %>% filter(SYMBOL %in% symbol.dict$Input)


symbols_to_add<-gene_symbols[which(!gene_symbols$SYMBOL %in% c(gene_ids$gene_symbol, gene_info$gene_symbol)),]

for(gene_id_update in na.omit(intersect(symbols_to_add$ENTREZID, gene_ids$ENTREZID))){
  gene_ids$gene_symbol[which(gene_ids$ENTREZID==gene_id_update)]<-symbols_to_add$SYMBOL[which(symbols_to_add$ENTREZID==gene_id_update)]
}
symbols_to_add<-symbols_to_add %>% filter(!ENTREZID %in% intersect(symbols_to_add$ENTREZID, gene_ids$ENTREZID))

for(gene_id_update in na.omit(intersect(symbols_to_add$ENTREZID, gene_info$ENTREZID))){
  gene_info$gene_symbol[which(gene_info$ENTREZID==gene_id_update)]<-symbols_to_add$SYMBOL[which(symbols_to_add$ENTREZID==gene_id_update)]
}
symbols_to_add<-symbols_to_add %>% filter(!ENTREZID %in% na.omit(intersect(symbols_to_add$ENTREZID, gene_info$ENTREZID)))

gene_symbols<-data.frame(gene_id=NA,
                         gene_symbol=symbols_to_add$SYMBOL,
                         ENTREZID=symbols_to_add$ENTREZID)

gene_dictionary<-rbind(gene_info, gene_ids, gene_symbols)

rm(gene_ids);rm(gene_info);rm(gene_symbols);rm(gene_symbols.entrez);rm(gene_symbols.tcga)
rm(symbol.dict);rm(symbol.info);rm(symbol.to.id);rm(symbols_to_add)
rm(whole.gene.entrez);rm(whole.gene.symbol);
rm(gene_id_update);rm(na.ind);rm(protein_coding_gene.ind);rm(symbol.update);rm(tumor.ind);rm(adc.ind)
rm(mad.exclusion);rm(qnorm);rm(series);rm(i)

dataset.name<-c("TCGA-LUAD", "GSE11969", "GSE13213", "GSE26939", "GSE31210", "GSE72094")

### Total data lists
expression.data.list<-list(expression_TCGA, expression_GSE11969, expression_GSE13213, 
                           expression_GSE26939, expression_GSE31210, expression_GSE72094)

rm(expression_TCGA);rm(expression_GSE11969); rm(expression_GSE13213);
rm(expression_GSE26939); rm(expression_GSE31210); rm(expression_GSE72094)
rm(temp); rm(labels)

for(i in c(2,3)){
  expression.data.list[[i]]<-
    expression.data.list[[i]][, -match(setdiff(colnames(expression.data.list[[i]]),
                                               gene_dictionary$gene_symbol),
                                       colnames(expression.data.list[[i]]))]
}


### negative data normalization for the datasets from Agilent platform
which(sapply(expression.data.list, min)<0) # Agilent platform: 2,3
expression.data.list[[2]]<-t(apply(expression.data.list[[2]], MARGIN=1, FUN=function(x){
  return((x- min(x)) /(max(x)-min(x)))
})*15)
expression.data.list[[3]]<-t(apply(expression.data.list[[3]], MARGIN=1, FUN=function(x){
  return((x- min(x)) /(max(x)-min(x)))
})*15)



clinical.annotation.list<-list(clinical_TCGA, clinical_GSE11969, clinical_GSE13213,
                               clinical_GSE26939, clinical_GSE31210, clinical_GSE72094)

names(clinical.annotation.list)<-names(expression.data.list)<-dataset.name

rm(clindata);rm(clinical);
rm(clinical_TCGA); rm(clinical_GSE11969); rm(clinical_GSE13213);
rm(clinical_GSE26939); rm(clinical_GSE31210); rm(clinical_GSE72094)
                                                                        
# datasets containing particular genes
clusterExport(clu, "gene_dictionary")
clusterExport(clu, "clinical.annotation.list")
saveRDS(expression.data.list, "expression.data.list.RDS")
clusterEvalQ(clu, {
  expression.data.list<-readRDS("expression.data.list.RDS")
  NULL
})
#clusterExport(clu, "expression.data.list")
gene_dictionary$dataset_info<-pbsapply(c(1:dim(gene_dictionary)[1]), FUN=function(k){
  gene_containing_dataset_ind<-which(sapply(c(1:length(expression.data.list)), FUN=function(i){
    if(length(na.omit(match(gene_dictionary[k,], colnames(expression.data.list[[i]]))))!=0){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }))
  return(paste0(gene_containing_dataset_ind, collapse="|"))
}, cl=clu)
