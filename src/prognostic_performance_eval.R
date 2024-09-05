score.type="risk score"
Dataset<-Marker_set<-C.index<-P.value<-n.patients<-c()
for(gset.ind in 1:length(test.markers)){
  # gene set 선정
  geneset<-test.markers[[gset.ind]]
  
  for(dset.ind in 1:length(evidence.datasets)){
    # leave-one-out 할 데이터셋 선정
    dataset.leaveoneout<-evidence.datasets[dset.ind]
    
    # leave-one-dataset-out 기반 summarized gene의 univariate Cox HR 도출
    summary_surv_res_geneset<-do.call(rbind, pblapply(match(geneset, gene_dictionary$gene_symbol), FUN=survival.marker.stat.summary, 
                                                      discovery.set=setdiff(evidence.datasets, evidence.datasets[dset.ind]), cl=clu))
    
    ## expression dataset 내 gene들의 index 도출
    dataset.ind<-match(dataset.leaveoneout, dataset.name)
    geneid.type.ind<-which(sapply(c(1:3), FUN=function(k){
      length(intersect(colnames(expression.data.list[[dataset.ind]]), gene_dictionary[,k]))
    })>1000)
    geneset.ind<-match(gene_dictionary[match(geneset, gene_dictionary$gene_symbol), geneid.type.ind], 
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
    
    if(score.type=="risk score"){
      # risk score 계산
      coefficients<-log(summary_surv_res_geneset$HR)
      names(coefficients)<-geneset
      
      ## 전체 환자에 대한 risk score 계산 후 condition 해당 patient index로 filter함
      score<-sapply(c(1:dim(expression.data.list[[dataset.ind]])[1]), FUN=function(k){
        patient.exp<-expression.data.list[[dataset.ind]][k,geneset.ind]
        score<-sum(na.omit(coefficients*patient.exp))
        return(score)
      })
      
    }
    
    if(score.type=="geometric mean difference score"){
      # geometric mean score 계산
      coefficients<-sign(log(summary_surv_res_geneset$HR))
      names(coefficients)<-geneset
      
      ## 전체 환자에 대한 geometric mean difference score 계산 후 condition 해당 patient index로 filter함
      score<-sapply(c(1:dim(expression.data.list[[dataset.ind]])[1]), FUN=function(k){
        patient.exp<-expression.data.list[[dataset.ind]][k,geneset.ind]
        score<-na.omit(coefficients*patient.exp)
        if(length(score)==0){
          return(0)
        }
        pos.score<-exp(mean(log(score[which(score>=0)])))
        neg.score<-exp(mean(log(-score[which(score<0)])))
        if(is.na(pos.score)){
          score<-(-neg.score)
        }else if(is.na(neg.score)){
          score<-pos.score
        }else{
          score<-pos.score-neg.score
        }
        
        return(score)
      })
      
    }
    
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
    
    Dataset<-c(Dataset, evidence.datasets[dset.ind])
    Marker_set<-c(Marker_set, names(test.markers)[[gset.ind]])
    C.index<-c(C.index, c.index)
    P.value<-c(P.value, km.pval)
    n.patients<-c(n.patients, length(patient.ind))
  }
  
}

df<-data.frame(Dataset, Marker_set, C.index, P.value, n.patients)

df$P.value<-sapply(c(1:dim(df)[1]), FUN=function(k){
  if(df$C.index[k]>=0.5){
    return(df$P.value[k])
  }else{
    return(1-df$P.value[k])
  }
})

cv.result<-do.call(rbind, lapply(unique(df$Marker_set), FUN=function(x){
  Marker_set=x
  mean=(df %>% filter(Marker_set==x))$C.index %>% mean
  sd=(df %>% filter(Marker_set==x))$C.index %>% sd
  return(data.frame(Marker_set, mean, sd))
}))

df.cv<-df

## 기존 마커 셋의 independent data에서의 성능 비교
pvals<-cind<-pvals.best<-hazard.ratio<-cox.pval<-c()
for(n in 1:length(test.markers)){
  
  target.gene.info<-summary_surv_res %>% filter(gene %in% test.markers[[n]])
  target.gene<-target.gene.info$gene
  
  ## 전체 환자에 대한 score 계산 후 condition 해당 patient index로 filter함
  dataset.ind<-match(validation.set, dataset.name)[1]
  geneid.type.ind<-which(sapply(c(1:3), FUN=function(k){
    length(intersect(colnames(expression.data.list[[dataset.ind]]), gene_dictionary[,k]))
  })>1000)
  geneset.ind<-match(gene_dictionary[match(target.gene, gene_dictionary$gene_symbol), geneid.type.ind], 
                     colnames(expression.data.list[[dataset.ind]]))
  
  if(score.type=="risk score"){
    
    # risk score 계산
    coefficients<-log(target.gene.info$HR)
    names(coefficients)<-target.gene
    score<-sapply(c(1:dim(expression.data.list[[dataset.ind]])[1]), FUN=function(k){
      patient.exp<-expression.data.list[[dataset.ind]][k,geneset.ind]
      score<-sum(na.omit(coefficients*patient.exp))
      return(score)
    })
    
  }else if(score.type=="geometric mean difference score"){
    # geometric mean score 계산
    coefficients<-sign(log(target.gene.info$HR))
    names(coefficients)<-target.gene.info$gene
    
    score<-sapply(c(1:dim(expression.data.list[[dataset.ind]])[1]), FUN=function(k){
      patient.exp<-expression.data.list[[dataset.ind]][k,geneset.ind]
      score<-na.omit(coefficients*patient.exp)
      if(length(score)==0){
        return(0)
      }
      pos.score<-exp(mean(log(score[which(score>=0)])))
      neg.score<-exp(mean(log(-score[which(score<0)])))
      if(is.na(pos.score)){
        score<-(-neg.score)
      }else if(is.na(neg.score)){
        score<-pos.score
      }else{
        score<-pos.score-neg.score
      }
      
      return(score)
    })
    
  }
  
  
  ## survival object generation
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
                 score=score[patient.ind])
  df<-df %>% filter(!time %in% na.terms) 
  df$event<-as.numeric(df$event)
  df$time<-as.numeric(df$time)
  df<-df %>% filter(time>=0) 
  
  surv.obj<-Surv(as.numeric(df$time), as.numeric(df$event))
  c.index<-1-concordance(surv.obj~df$score)$concordance
  cind<-c(cind, c.index)
  
  cox.model<-coxph(Surv(as.numeric(time), as.numeric(event))~., data=df) 
  HR<-exp(cox.model$coefficients)
  hazard.ratio<-c(hazard.ratio, HR)
  
  cox.pval<-c(cox.pval, summary(cox.model)$coefficients[1,5])
  
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
  pvals<-c(pvals, km.pval)
  
  
  #### best prognostic group stratification
  best.search<-list()
  count<-1
  for(threshold in sort(df$score)[ceiling(length(df$score)*0.25):floor(length(df$score)*0.75)]){
    df2=df
    df2$group=sapply(df$score, FUN=function(k){
      if(k>=threshold){
        return("high_score")
      }else{
        return("low_score")
      }
    })
    surv.obj<-Surv(as.numeric(df2$time), as.numeric(df2$event))
    fit.survival<-survfit(surv.obj~group, data=df2)
    km.pval<-surv_pvalue(fit.survival)$pval
    best.search[[count]]<-data.frame(threshold=threshold, p.value=km.pval)
    count<-count+1
  }
  best.search<-do.call(rbind, best.search)
  
  threshold<-best.search$threshold[which(best.search$p.value==min(best.search$p.value))]
  
  df$group=sapply(df$score, FUN=function(k){
    if(k>=median(best.search$threshold[which(best.search$p.value==min(best.search$p.value))])){
      return("high_score")
    }else{
      return("low_score")
    }
  })
  surv.obj<-Surv(as.numeric(df$time), as.numeric(df$event))
  fit.survival<-survfit(surv.obj~group, data=df)
  km.pval2<-surv_pvalue(fit.survival)$pval
  
  pvals.best<-c(pvals.best, km.pval2)
}

independent.res<-data.frame(geneset=names(test.markers), HR=hazard.ratio, c.index=cind, p.value=cox.pval, p.value.best=pvals.best)
