validation.set="TCGA-LUAD"


### identifying optimal number of marker
SYNBI_Entire=(whole.patients.significant.gene.info %>% arrange(-c.ind.min))$gene[1:81]
SYNBI_non_EGFR_Mut=(egfr.wt.significant.gene.info %>% arrange(-c.ind.min))$gene[1:83]
SYNBI_EGFR_Mut=(egfr.significant.gene.info %>% arrange(-c.ind.min))$gene[1:32]

Table2<-(egfr.significant.gene.info %>% arrange(-c.ind.min))[1:32,] %>%
  mutate(Coef=log(HR)) 

Table2$CI<-sapply(c(1:dim(Table2)[1]), FUN=function(i){
  paste0(c(round(Table2$lower[i],4), "-", round(Table2$upper[i], 4)), collapse="")
})

Table2<-data.frame(Gene=Table2$gene, Coefficient=round(Table2$Coef, 4),
                   HR=round(Table2$HR, 4), CI=Table2$CI, adj_p=round(Table2$adjusted.p, 4))

write.csv(Table2, file = "Table2.csv", row.names = F, col.names = T)

test.markers<-c(previous_markers, list(SYNBI_Entire=SYNBI_Entire,
                                       SYNBI_non_EGFR_Mut=SYNBI_non_EGFR_Mut,
                                       SYNBI_EGFR_Mut=SYNBI_EGFR_Mut))


#### EGFR mutant
condition<-"egfr_mut:EGFR-Mut"
each.gene.survival.res<-survival.marker.stat.individual(input_query=condition)
saveRDS(each.gene.survival.res, "each.gene.survival.res.RDS")

clusterEvalQ(clu, {
  each.gene.survival.res<-readRDS('each.gene.survival.res.RDS')
  NULL
})
file.remove('each.gene.survival.res.RDS')
summary_surv_res=summary_surv_res.egfr_mutant

source("prognostic_performance_eval.R")


df.cv.egfr_mut<-df.cv
cv.result.egfr_mut<-cv.result
independent.validation.egfr.mut<-independent.res

Figure4a<-ggplot(df.cv.egfr_mut, aes(x=factor(Marker_set, levels = unique(Marker_set)),
               y=C.index))+
  geom_bar(data=cv.result.egfr_mut, aes(x=factor(Marker_set, levels=unique(Marker_set)),
                               y=mean), stat='identity', fill="gray")+
  geom_boxplot(width=0.4, outlier.size=0.05)+
  geom_point(size=1)+
  geom_hline(yintercept=0.5, linetype="dashed")+
  scale_y_continuous(limits = c(0, .95))+
  geom_text(data=cv.result.egfr_mut, aes(y=0.95, label=round(mean, 3)), position=position_dodge(width=.75), size=2)+
  theme_bw()+
  scale_x_discrete(labels=sapply(unique(df.cv.egfr_mut$Marker_set), FUN=function(x){
    if(x=="SYNBI_EGFR_Mut"){
      return(expression(bold("SYNBI_EGFR_Mut")))
    }else
      return(x)
  }))+
  labs(x="Marker set", y="C-index")+
  ggtitle('EGFR-mutated patients, cross-validation')+
  coord_flip()+
  theme(axis.text=element_text(size=6), plot.title = element_text(size=8), axis.title=element_text(size=8))

tiff(file="Figure4a.tiff", width = 13, height = 7.5, units = "cm", res=300)
Figure4a
dev.off()


Figure5a<-ggplot(independent.validation.egfr.mut, 
       aes(x=factor(geneset,
                    levels=(independent.validation.egfr.mut %>% arrange(c.index))$geneset),
           y=c.index))+
  geom_bar(stat='identity', fill="gray")+
  xlab("")+ylab("C-index")+
  scale_y_continuous(limits = c(0, 0.8))+
  geom_text(data=independent.validation.egfr.mut, aes(y=c.index, label=round(c.index, 3)), 
            position=position_dodge(width=.75), size=1.8)+
  scale_x_discrete(labels=sapply((independent.validation.egfr.mut %>% arrange(c.index))$geneset, FUN=function(x){
    if(x=="SYNBI_EGFR_Mut"){
      return(expression(bold("SYNBI_EGFR_Mut")))
    }else
      return(x)
  }))+
  coord_flip()+
  theme_bw()+
  ggtitle("EGFR-mutated patients, independent validation \n(combined risk coefficients)")+
  theme(axis.text=element_text(size=6), plot.title = element_text(size=8), axis.title=element_text(size=8))

tiff(file="Figure5a.tiff", width = 10, height = 7.5, units = "cm", res=300)
Figure5a
dev.off()

df.surv.independent<-df
df.surv.independent$group=sapply(df.surv.independent$score, FUN=function(k){
  if(k>median(df.surv.independent$score)){
    return("high_score")
  }else{
    return("low_score")
  }
})
surv.obj<-Surv(12*as.numeric(df.surv.independent$time)/365, as.numeric(df.surv.independent$event))
fit.survival<-survfit(surv.obj~group, data=df.surv.independent)
Figure5c<-ggsurvplot(fit.survival, pval = T, conf.int = T, risk.table = "nrisk_cumcensor")+xlab("Time (Month)")

tiff(file="Figure5c.tiff", width=18, height=15, units="cm", res=300)
Figure5c+
  xlab("Time (Months)")
dev.off()

df.surv.independent<-list(df.surv.independent)

#### Whole patients
condition=""
each.gene.survival.res<-survival.marker.stat.individual(input_query=condition)
saveRDS(each.gene.survival.res, "each.gene.survival.res.RDS")

clusterEvalQ(clu, {
  each.gene.survival.res<-readRDS('each.gene.survival.res.RDS')
  NULL
})
file.remove('each.gene.survival.res.RDS')
summary_surv_res=summary_surv_res.whole_patients
source("prognostic_performance_eval.R")

df.cv.entire_patient<-df.cv
cv.result.entire.patient<-cv.result


independent.validation.entire.patients<-independent.res

Figure4c<-ggplot(df.cv.entire_patient, aes(x=factor(Marker_set, levels = unique(Marker_set)),
                  y=C.index))+
  geom_bar(data=cv.result.entire.patient, aes(x=factor(Marker_set, levels=unique(Marker_set)),
                               y=mean), stat='identity', fill="gray")+
  geom_boxplot(width=0.4, outlier.size=0.05)+
  geom_point(size=1)+
  geom_hline(yintercept=0.5, linetype="dashed")+
  scale_y_continuous(limits = c(0, 0.95)) +
  geom_text(data=cv.result.entire.patient, aes(y=0.95, label=round(mean, 3)), position=position_dodge(width=.75), size=2)+
  theme_bw()+
  labs(x="Marker set", y="C-index")+
  coord_flip()+
  ggtitle('Entire patients, cross-validation')+
  theme(axis.text=element_text(size=6), plot.title = element_text(size=8), axis.title=element_text(size=8))

tiff(file="Figure4c.tiff", width = 13, height = 7.5, units = "cm", res=300)
Figure4c
dev.off()

##  EGFR wild type
condition<-"egfr_mut:EGFR-WT"
each.gene.survival.res<-survival.marker.stat.individual(input_query=condition)
saveRDS(each.gene.survival.res, "each.gene.survival.res.RDS")

clusterEvalQ(clu, {
  each.gene.survival.res<-readRDS('each.gene.survival.res.RDS')
  NULL
})
file.remove('each.gene.survival.res.RDS')
summary_surv_res=summary_surv_res.egfr_wt
source("prognostic_performance_eval.R")

df.cv.egfr_wt<-df.cv
cv.result.egfr_wt<-cv.result
independent.validation.egfr_wt<-independent.res


Figure4b<-ggplot(df.cv.egfr_wt, aes(x=factor(Marker_set, levels = unique(Marker_set)),
                  y=C.index))+
  geom_bar(data=cv.result.egfr_wt, aes(x=factor(Marker_set, levels=unique(Marker_set)),
                               y=mean), stat='identity', fill="gray")+
  geom_boxplot(width=0.4, outlier.size = 0.05)+
  geom_point(size=1)+
  geom_hline(yintercept=0.5, linetype="dashed")+
  scale_y_continuous(limits = c(0, 0.95)) +
  geom_text(data=cv.result.egfr_wt, aes(y=0.95, label=round(mean, 3)), position=position_dodge(width=.75), size=2)+
  theme_bw()+
  labs(x="Marker set", y="C-index")+
  coord_flip()+
  ggtitle('Non-EGFR-mutated patients, cross-validation')+
  theme(axis.text=element_text(size=6), plot.title = element_text(size=8), axis.title=element_text(size=8))

tiff(file="Figure4b.tiff", width = 13, height = 7.5, units = "cm", res=300)
Figure4b
dev.off()

### Performance evaluation using the markers with available risk coefficients
previous_markers_coef<-read.csv("previous_marker_with_coef.CSV")
previous_markers_coef<-lapply(unique(previous_markers_coef$Study), FUN=function(x){
  previous_markers_coef %>% filter(Study==x)
})

names(previous_markers_coef)<-sapply(previous_markers_coef, FUN=function(x){x$Study[1]})

evidence.datasets<-c("TCGA-LUAD", evidence.datasets)

Dataset<-Marker_set<-C.index<-P.value<-c()
for(gset.ind in 1:length(previous_markers_coef)){
  
  for(dset.ind in 1:length(evidence.datasets)){# test 할 데이터셋
    ## expression dataset 내 gene들의 index 도출
    dataset.ind<-match(evidence.datasets[dset.ind], dataset.name)
    
    geneid.type.ind<-which(sapply(c(1:3), FUN=function(k){
      length(intersect(colnames(expression.data.list[[dataset.ind]]), gene_dictionary[,k]))
    })>1000)
    geneset.ind<-match(gene_dictionary[match(previous_markers_coef[[gset.ind]]$Gene, 
                                             gene_dictionary$gene_symbol), geneid.type.ind], 
                       colnames(expression.data.list[[dataset.ind]]))
    
    
    # condition 해당하는 patient index filtering
    clinical.data<-clinical.annotation.list[[dataset.ind]]
    if(condition!=""){
      category<-sapply(unlist(strsplit((condition %>% strsplit("+", fixed=T) %>% unlist)[1], "&", fixed=T)),
                       FUN=function(x){
                         unlist(strsplit(x, ":", fixed=T))[1]
                       })
      confounding.variables<-(condition %>% strsplit("+", fixed=T) %>% unlist)[2] %>% strsplit("&", fixed=T) %>% unlist
      
      
      if(!(length(category)==0||sum(is.na(category))>0)){
        label<-lapply(unlist(strsplit((condition %>% strsplit("+", fixed=T) %>% unlist)[1], "&")), FUN=function(x){
          strsplit(unlist(strsplit(x, ":"))[2], "|", fixed=T) %>% unlist
        })
      }
      
      ## confounding variable의 경우, confounding variable의 label들이 2개 이상 존재하는 것(valid variables)들만 활용함.
      if(!(length(confounding.variables)==0 || is.na(confounding.variables))){
        confounding.variable.validity<-sapply(confounding.variables, FUN=function(x){ 
          # if confounding variable is contained by category to analyze, it should contain 2 or more valid values to build Cox PH model
          if(is.na(match(x, category))){
            return(TRUE)
          }else{
            return(length(label.to.analyze[[match(x, category)]])>1)
          }
        })
        
        confounding.variables<-confounding.variables[which(confounding.variable.validity)]
      }
      
      if(length(confounding.variables)==0){
        confounding.variables=NA
      }
      
      patient.ind.each.condition<-lapply(c(1:length(category)), FUN=function(n){
        which(clinical.data[,match(category[n], colnames(clinical.data))] %in% label[[n]])
      })
      
      if(!is.na(confounding.variables)){
        patient.ind.valid.confounding<-lapply(c(1:length(confounding.variables)), FUN=function(n){
          which(!clinical.data[,match(confounding.variables[n], colnames(clinical.data))] %in% na.terms)
        })
        patient.ind<-Reduce(intersect, c(patient.ind.each.condition, patient.ind.valid.confounding))
      }else{
        patient.ind<-Reduce(intersect, patient.ind.each.condition)
      }
      
      
    }else{
      patient.ind<-c(1:dim(clinical.data)[1])
    }
    
    if(length(patient.ind)<min.patient.num){
      next()
    }
    
    # risk score 계산
    coefficients<-previous_markers_coef[[gset.ind]]$Coef
    names(coefficients)<-previous_markers_coef[[gset.ind]]$Gene
    score<-sapply(c(1:dim(expression.data.list[[dataset.ind]])[1]), FUN=function(k){
      patient.exp<-expression.data.list[[dataset.ind]][k,geneset.ind]
      score<-sum(na.omit(coefficients*patient.exp))
      return(score)
    })
    
    if(sum(score)==0){
      next()
    }
    
    ## survival object generation
    df<-data.frame(event=ifelse(clinical.data$vital_status[patient.ind]=="Dead", 1, 0),
                   time=clinical.data$survival_time[patient.ind],
                   score=score[patient.ind])
    df<-df %>% filter(!time %in% na.terms) 
    df$event<-as.numeric(df$event)
    df$time<-as.numeric(df$time)
    df<-df %>% filter(time>=0) 
    
    surv.obj<-Surv(as.numeric(df$time), as.numeric(df$event))
    c.index<-1-concordance(surv.obj~df$score)$concordance
    
    df$group=sapply(df$score, FUN=function(k){
      if(k>=median(df$score)){
        return("high_score")
      }else{
        return("low_score")
      }
    })
    surv.obj<-Surv(as.numeric(df$time), as.numeric(df$event))
    fit.survival<-survfit(surv.obj~group, data=df)
    km.pval<-surv_pvalue(fit.survival)$pval
    
    #### independent dataset
    if(dset.ind==1){
      df.surv.independent<-c(df.surv.independent, list(df))
    }
    
    if(gset.ind==length(test.markers)){
      surv.plot[[count]]<-ggsurvplot(fit.survival, data=df, pval=TRUE, risk.table = TRUE,  risk.table.y.text.col=FALSE,
                                     title=evidence.datasets[dset.ind], ggtheme=theme_classic2(base_size=10))
      
      count<-count+1
    }
    
    Dataset<-c(Dataset, evidence.datasets[dset.ind])
    Marker_set<-c(Marker_set, names(previous_markers_coef)[gset.ind])
    C.index<-c(C.index, c.index)
    P.value<-c(P.value, km.pval)
  }
  
}

evidence.datasets<-evidence.datasets[-1]

df.performance.with.prev.coef<-data.frame(Dataset, Marker_set, C.index, P.value)
df.performance.with.prev.coef.discovery<-df.performance.with.prev.coef %>% filter(Dataset!="TCGA-LUAD")
cv.result.prev.coef.discovery<-do.call(rbind, lapply(unique(df.performance.with.prev.coef.discovery$Marker_set), FUN=function(x){
  Marker_set=x
  mean=(df.performance.with.prev.coef.discovery %>% filter(Marker_set==x))$C.index %>% mean
  sd=(df.performance.with.prev.coef.discovery %>% filter(Marker_set==x))$C.index %>% sd
  return(data.frame(Marker_set, mean, sd))
}))

FigureS2<-ggplot(df.performance.with.prev.coef.discovery, 
                  aes(x=factor(Marker_set, levels = unique(Marker_set)),
                    y=C.index))+
  geom_bar(data=cv.result.prev.coef.discovery, aes(x=factor(Marker_set, levels=unique(Marker_set)),
                               y=mean), stat='identity', fill="lightgray")+
  geom_text(data=cv.result.prev.coef.discovery, aes(y=0.9, label=round(mean, 3)), position=position_dodge(width=.75), size=4)+
  geom_boxplot(width=0.4)+
  geom_point()+
  geom_hline(yintercept=0.5, linetype="dashed")+
  theme_classic()+
  scale_y_continuous(limits = c(0, 0.9)) +
  labs(x="Marker set", y="Concordance index")+
  coord_flip()+
  ggtitle("EGFR-mutated patients, Discovery sets \n(previously given risk coefficients)")+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=10))

tiff(file="FigureS2.tiff", width = 15, height = 10, units = "cm", res=300)
FigureS2
dev.off()

df.performance.with.prev.coef.validation<-df.performance.with.prev.coef %>% filter(Dataset=="TCGA-LUAD")
df.performance.with.prev.coef.validation<-
  rbind(df.performance.with.prev.coef.validation,
        data.frame(Dataset="TCGA-LUAD", Marker_set="SYNBI_EGFR_Mut",
                   C.index=independent.validation.egfr.mut$c.index[28],
                   P.value=independent.validation.egfr.mut$p.value[28]))

Figure5b<-ggplot(df.performance.with.prev.coef.validation, 
                 aes(x=factor(Marker_set,
                              levels=(df.performance.with.prev.coef.validation %>% arrange(C.index))$Marker_set),
                     y=C.index))+
  geom_bar(stat='identity', fill="gray")+
  xlab("Marker set")+ylab("C-index")+
  scale_y_continuous(limits = c(0, 0.8))+
  geom_text(data=df.performance.with.prev.coef.validation, aes(y=C.index, label=round(C.index, 3)), 
            position=position_dodge(width=.75), size=1.8)+
  scale_x_discrete(labels=sapply((df.performance.with.prev.coef.validation %>% arrange(C.index))$Marker_set, FUN=function(x){
    if(x=="SYNBI_EGFR_Mut"){
      return(expression(bold("SYNBI_EGFR_Mut")))
    }else
      return(x)
  }))+
  coord_flip()+
  theme_bw()+
  ggtitle("EGFR-mutated patients, independent validation \n(previously given risk coefficients)")+
  theme(axis.text=element_text(size=6), plot.title = element_text(size=8), axis.title=element_text(size=8))

tiff(file="Figure5b.tiff", width = 10, height = 5, units = "cm", res=300)
Figure5b
dev.off()


# Multivariate Cox regression with confounding variables: age, gender, and pathological stages.
condition<-"egfr_mut:EGFR-Mut+age&gender&ajcc_pathologic_stage"
target.gene.info<-summary_surv_res %>% filter(gene %in% test.markers$SYNBI_EGFR_Mut)
target.gene<-target.gene.info$gene


cox.adjusted<-lapply(egfr.datasets, FUN=function(x){
  dataset.ind<-match(x, dataset.name)
  geneid.type.ind<-which(sapply(c(1:3), FUN=function(k){
    length(intersect(colnames(expression.data.list[[dataset.ind]]), gene_dictionary[,k]))
  })>1000)
  geneset.ind<-match(gene_dictionary[match(target.gene, gene_dictionary$gene_symbol), geneid.type.ind], 
                     colnames(expression.data.list[[dataset.ind]]))
  
  
  # risk score calculation
  coefficients<-log(target.gene.info$HR)
  names(coefficients)<-target.gene
  score<-sapply(c(1:dim(expression.data.list[[dataset.ind]])[1]), FUN=function(k){
    patient.exp<-expression.data.list[[dataset.ind]][k,geneset.ind]
    score<-sum(na.omit(coefficients*patient.exp))
    return(score)
  })
  
  ## survival object generation
  # conditional patient index filtering
  clinical.data<-clinical.annotation.list[[dataset.ind]]
  if(condition!=""){
    category<-sapply(unlist(strsplit((condition %>% strsplit("+", fixed=T) %>% unlist)[1], "&", fixed=T)),
                     FUN=function(x){
                       unlist(strsplit(x, ":", fixed=T))[1]
                     })
    confounding.variables<-(condition %>% strsplit("+", fixed=T) %>% unlist)[2] %>% strsplit("&", fixed=T) %>% unlist
    
    
    if(!(length(category)==0||sum(is.na(category))>0)){
      label<-lapply(unlist(strsplit((condition %>% strsplit("+", fixed=T) %>% unlist)[1], "&")), FUN=function(x){
        strsplit(unlist(strsplit(x, ":"))[2], "|", fixed=T) %>% unlist
      })
    }
    
    ## confounding variable의 경우, confounding variable의 label들이 2개 이상 존재하는 것(valid variables)들만 활용함.
    if(!(length(confounding.variables)==0 || is.na(confounding.variables))){
      confounding.variable.validity<-sapply(confounding.variables, FUN=function(x){ 
        # if confounding variable is contained by category to analyze, it should contain 2 or more valid values to build Cox PH model
        if(is.na(match(x, category))){
          return(TRUE)
        }else{
          return(length(label.to.analyze[[match(x, category)]])>1)
        }
      })
      
      confounding.variables<-confounding.variables[which(confounding.variable.validity)]
    }
    
    if(length(confounding.variables)==0){
      confounding.variables=NA
    }
    
    patient.ind.each.condition<-lapply(c(1:length(category)), FUN=function(n){
      which(clinical.data[,match(category[n], colnames(clinical.data))] %in% label[[n]])
    })
    
    if(!is.na(confounding.variables[1])){
      patient.ind.valid.confounding<-lapply(c(1:length(confounding.variables)), FUN=function(n){
        which(!clinical.data[,match(confounding.variables[n], colnames(clinical.data))] %in% na.terms)
      })
      patient.ind<-Reduce(intersect, c(patient.ind.each.condition, patient.ind.valid.confounding))
      
    }else{
      patient.ind<-Reduce(intersect, patient.ind.each.condition)
    }
    
  }else{
    patient.ind<-c(1:dim(clinical.data)[1])
  }
  
  df<-data.frame(event=ifelse(clinical.data$vital_status[patient.ind]=="Dead", 1, 0),
                 time=clinical.data$survival_time[patient.ind],
                 score=score[patient.ind],
                 clinical.data[patient.ind,confounding.variables])
  
  df$ajcc_pathologic_stage<-sapply(df$ajcc_pathologic_stage, FUN=function(x){
    if(x %in% c("Stage I", "Stage IA", "Stage IB")){
      return("Stage I")
    }else if(x %in% c("Stage IIA", "Stage IIB", "Stage II")){
      return("Stage II")
    }else if(x %in% c("Stage IIIA", "Stage IIIB", "Stage III")){
      return("Stage III")
    }else if(x %in% "Stage IV"){
      return("Stage IV")
    }else{
      return(x)
    }
  })
  df<-df %>% filter(!time %in% na.terms) 
  df$event<-as.numeric(df$event)
  df$time<-as.numeric(df$time)
  df$age<-as.numeric(df$age)
  df<-df %>% filter(time>=0) 
  
  cox.model<-coxph(Surv(as.numeric(time), as.numeric(event))~., data=df) 
  return(cox.model)
})

