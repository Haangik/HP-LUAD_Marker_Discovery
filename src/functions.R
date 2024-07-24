survival.marker.stat.individual<-function(input_query, condition.disjoint=F, egfr.shuffle=F, other.shuffle=F, shuffle.portion=0.7){
  category.to.analyze<-sapply(unlist(strsplit((input_query %>% strsplit("+", fixed=T) %>% unlist)[1], "&", fixed=T)),
                              FUN=function(x){
                                unlist(strsplit(x, ":", fixed=T))[1]
                              })
  
  confounding.variables<-(input_query %>% strsplit("+", fixed=T) %>% unlist)[2] %>% strsplit("&", fixed=T) %>% unlist
  
  if(!(length(category.to.analyze)==0||sum(is.na(category.to.analyze))>0)){
    label.to.analyze<-lapply(unlist(strsplit((input_query %>% strsplit("+", fixed=T) %>% unlist)[1], "&")), FUN=function(x){
      strsplit(unlist(strsplit(x, ":"))[2], "|", fixed=T) %>% unlist
    })
  }
  
  if(!(length(confounding.variables)==0 || is.na(confounding.variables))){
    confounding.variable.validity<-sapply(confounding.variables, FUN=function(x){ 
      # if confounding variable is contained by category to analyze, it should contain 2 or more valid values to build Cox PH model
      if(is.na(match(x, category.to.analyze))){
        return(TRUE)
      }else{
        return(length(label.to.analyze[[match(x, category.to.analyze)]])>1)
      }
    })
    
    confounding.variables<-confounding.variables[which(confounding.variable.validity)]
  }
  if(length(confounding.variables)==0){
    confounding.variables=NA
  }
  
  
  cat(paste0("Target patient condition information \n",
             ifelse((is.na(category.to.analyze)||length(category.to.analyze)==0),
                    "No specified clinical condition", 
                    paste0(sapply(c(1:length(category.to.analyze)), FUN=function(k){
                      paste0(category.to.analyze[k], ": ", paste0(label.to.analyze[[k]], collapse = ", "))
                    }), collapse="\n")),
             "\n",
             "Confounding variable: ",
             ifelse(is.na(confounding.variables), "None", paste0(confounding.variables, collapse=", ")),
             "\n"
  ))
  
  dataset.patient.ind<-lapply(c(1:length(dataset.name)), FUN=function(k){
    if(length(category.to.analyze)==0||is.na(category.to.analyze)){
      each.condition.patient.ind<-list()
      each.condition.patient.ind[[1]]<-c(1:dim(clinical.annotation.list[[k]])[1])
    }else{
      each.condition.patient.ind<-lapply(c(1:length(category.to.analyze)), FUN=function(i){
        if(is.na(match(category.to.analyze[i], colnames(clinical.annotation.list[[k]])))){
          patient.ind<-0
        }else{
          labels<-clinical.annotation.list[[k]][, match(category.to.analyze[i], colnames(clinical.annotation.list[[k]]))]
          patient.ind<-which(labels %in% label.to.analyze[[i]])
        }
        return(patient.ind)
      })
    }
    
    if(!is.na(confounding.variables)){
      patient.with.confounding.ind<-lapply(confounding.variables, FUN=function(x){
        if(is.na(match(x, colnames(clinical.annotation.list[[k]])))){
          patient.ind<-0
        }else{
          labels<-clinical.annotation.list[[k]][,match(x, colnames(clinical.annotation.list[[k]]))]
          patient.ind<-which(!labels %in% na.terms)  
        }
        return(patient.ind)
      })
      each.condition.patient.ind<-c(each.condition.patient.ind, patient.with.confounding.ind)
    }
    
    common.patient.ind<-Reduce(intersect, each.condition.patient.ind)
    
    if(length(common.patient.ind)<min.patient.num){
      common.patient.ind<-integer(0)
      return(common.patient.ind)
    }
    
    vital.status<-clinical.annotation.list[[k]]$vital_status[common.patient.ind]
    if(min(length(which(vital.status=="Alive")), length(which(vital.status=="Dead")))<min.other.status){
      common.patient.ind<-integer(0)
      return(common.patient.ind)
    }
    return(common.patient.ind)
  })
  
  ### effects of patients random sampling
  if(egfr.shuffle==T){
    dataset.patient.ind<-lapply(c(1:length(dataset.patient.ind)), FUN=function(k){
      num.patient=ceiling(length(dataset.patient.ind[[k]])*shuffle.portion)
      return(sort(sample(dataset.patient.ind[[k]], num.patient)))
    })
  }
  
  if(other.shuffle==T){
    if(condition.disjoint==T){
      # non-mutant sampling
      dataset.patient.ind<-lapply(c(1:length(dataset.patient.ind)), FUN=function(k){
        num.patient=ceiling(length(dataset.patient.ind[[k]])*shuffle.portion)
        return(sort(sample(setdiff(c(1:dim(expression.data.list[[k]])[1]), dataset.patient.ind[[k]]),
                           num.patient, replace=T)))
      })
    }else{
      # entire patient sampling
      dataset.patient.ind<-lapply(c(1:length(dataset.patient.ind)), FUN=function(k){
        num.patient=ceiling(length(dataset.patient.ind[[k]])*shuffle.portion)
        return(sort(sample(dim(expression.data.list[[k]])[1], num.patient)))
      })
    }
  }
  
  
  cat(paste0("Total available patients: ", 
             sum(sapply(dataset.patient.ind, length)),
             " patients from ",
             length(which(sapply(dataset.patient.ind, length)!=0)),
             " dataset(s)\n"
  ))
  
  
  
  clusterExport(clu, c("dataset.patient.ind", "confounding.variables"), envir = environment())
  # Univariate Cox regression for each gene and dataset
  individual_dataset_surv_res<-pblapply(c(1:dim(gene_dictionary)[1]), FUN=function(g.ind){
    ds.ind<-intersect(unlist(strsplit(gene_dictionary$dataset_info[g.ind], "|", fixed=T)), 
                      which(sapply(dataset.patient.ind, length)>0)) %>% as.numeric
    
    if(length(ds.ind)==0){
      surv.res<-data.frame(dataset=character(), HR=numeric(),
                           lower=numeric(), upper=numeric(), p.value=numeric(), c.ind=numeric())
    }else{
      surv.res<-do.call(rbind, lapply(ds.ind, FUN=function(ds.ind){
        cat(ds.ind)
        patient.ind<-dataset.patient.ind[[ds.ind]]
        
        if(!is.na(confounding.variables)){
          surv.df<-data.frame(time=clinical.annotation.list[[ds.ind]]$survival_time[patient.ind],
                              event=ifelse(clinical.annotation.list[[ds.ind]]$vital_status[patient.ind]=="Dead", 1, 0),
                              expr=expression.data.list[[ds.ind]][patient.ind, 
                                                                  na.omit(match(gene_dictionary[g.ind,],
                                                                                colnames(expression.data.list[[ds.ind]])))],
                              clinical.annotation.list[[ds.ind]][patient.ind,
                                                                 match(confounding.variables,
                                                                       colnames(clinical.annotation.list[[ds.ind]]))]) %>% 
            filter(!time %in% na.terms) %>% filter(!event %in% na.terms) 
          colnames(surv.df)<-c('time', 'event', 'expr', confounding.variables)
        }else{
          surv.df<-data.frame(time=clinical.annotation.list[[ds.ind]]$survival_time[patient.ind],
                              event=ifelse(clinical.annotation.list[[ds.ind]]$vital_status[patient.ind]=="Dead", 1, 0),
                              expr=expression.data.list[[ds.ind]][patient.ind, 
                                                                  na.omit(match(gene_dictionary[g.ind,],
                                                                                colnames(expression.data.list[[ds.ind]])))]) %>% 
            filter(!time %in% na.terms) %>% filter(!event %in% na.terms)
        }
        
        surv.obj<-Surv(as.numeric(surv.df$time), as.numeric(surv.df$event))
        tryCatch({
          cox_ph_model<-coxph(Surv(as.numeric(time), as.numeric(event))~., data=surv.df)  
        },
        error=function(e){})
        if(!exists("cox_ph_model")){
          res<-data.frame(dataset=character(), HR=numeric(), lower=numeric(), upper=numeric(), p.value=numeric(), c.ind=numeric())
        }else{
          model.summary<-summary(cox_ph_model)
          HR<-model.summary$conf.int[1, 1]
          lower<-model.summary$conf.int[1, 3]
          upper<-model.summary$conf.int[1, 4]
          p.value<-model.summary$coefficients[1, 5]
          c.ind=model.summary$concordance[1]
          res<-data.frame(dataset=dataset.name[ds.ind], HR, lower, upper, p.value, c.ind)
        }
        
        return(res)
      }))
    }
    return(surv.res)
  }, cl=clu)
  
  
  return(individual_dataset_surv_res)
}


survival.marker.stat.summary<-function(k, discovery.set, genes=NULL){
  if(is.null(genes)){
    ## gene name
    if(is.na(gene_dictionary$gene_symbol[k])){
      gene<-gene_dictionary$ENTREZID[k]
    }else{
      gene<-gene_dictionary$gene_symbol[k]
    }
    
    if(gene==""){
      gene<-gene_dictionary$gene_id[k]
    }
  }else{
    gene<-genes[k]
  }
  
  surv.res<-each.gene.survival.res[[k]] %>% filter(dataset %in% discovery.set)
  if(is.null(surv.res) || dim(surv.res)[1]==0){
    summary.res<-data.frame(gene=gene,
                            HR=NA,
                            min.abs.log.HR=NA,
                            lower=NA,
                            upper=NA,
                            p.value=NA,
                            p.value.max=NA,
                            n.evidence=0,
                            c.ind.min=NA,
                            c.ind.mean=NA,
                            pval.Q=NA,
                            evidence_set=NA)
  }else{
    m<-metagen(HR=log(surv.res$HR),
               lower= log(surv.res$lower),
               upper= log(surv.res$upper),
               sm="HR", fixed=F, random=T,
               studlab = rownames(surv.res),
               method.tau="DL",
               method.random.ci = "classic")
    summary.res<-summary(m)
    summary.res<-data.frame(gene=gene,
                            HR=exp(summary.res$random$TE),
                            min.abs.log.HR=min(abs(log2(surv.res$HR))),
                            lower=exp(summary.res$random$lower),
                            upper=exp(summary.res$random$upper),
                            p.value=summary.res$random$p,
                            p.value.max=max(surv.res$p.value),
                            n.evidence=length(surv.res$dataset),
                            c.ind.min=min(surv.res$c.ind),
                            c.ind.mean=mean(surv.res$c.ind),
                            pval.Q=summary.res$pval.Q,
                            evidence_set=paste0(surv.res$dataset, collapse = "|"))
  }
  
  return(summary.res)
}

inhouse.annotation<-function(n, target.gene){
  gene.in.mod<-f2g$Entrez[which(f2g$MOD_ID==n)]
  gene.in.marker<-target.gene
  
  gene.in.both<-intersect(gene.in.mod, gene.in.marker)
  gene.in.mod.only<-setdiff(gene.in.mod, gene.in.both)
  gene.in.marker.only<-setdiff(gene.in.marker, gene.in.both)
  gene.in.none<-setdiff(total.genes, c(gene.in.mod, gene.in.marker))
  
  ftest.case<-fisher.test(matrix(c(length(gene.in.both), length(gene.in.mod.only),
                                   length(gene.in.marker.only), length(gene.in.none)), nrow=2))
  
  # Percentage of marker genes in module
  marker.in.mod=paste0(c(length(gene.in.both), "/", length(gene.in.both)+length(gene.in.marker.only), " (", 
                         round(length(gene.in.both)/(length(gene.in.both)+length(gene.in.marker.only)),3)*100, "%)"), collapse="")
 
  # percentage of module genes in marker
  mod.in.marker<-paste0(c(length(gene.in.both), "/", length(gene.in.mod), 
                          " (", round(length(gene.in.both)/length(gene.in.mod),3)*100, "%)"), collapse="")
  
  
  return(data.frame(mod.id=f$MOD_ID[n],
                    mod.db=f$DB[n],
                    mod.src=f$Source[n],
                    mod.name=f$Name[n],
                    annotation=f$Annotation[n],
                    pval=ftest.case$p.value,
                    OR=ftest.case$estimate,
                    marker.in.mod,
                    mod.in.marker
  ))
}

symbol.to.entrez<-function(symbol){
  gene_dictionary$ENTREZID[match(symbol, gene_dictionary$gene_symbol)]
}
entrez.to.symbol<-function(entrez){
  gene_dictionary$gene_symbol[match(entrez, gene_dictionary$ENTREZID)]
}
