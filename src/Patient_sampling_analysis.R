set.seed(123)
significant.marker.nonmut.sampling<-list()
median.favorable.nonmut.sampling<-median.unfavorable.nonmut.sampling<-c()
num.favorable.nonmut.sampling<-num.unfavorable.nonmut.sampling<-c()
c.index.nonmut.sampling<-c()
count<-1
while(TRUE){
  each.gene.survival.res<-survival.marker.stat.individual(input_query=condition, other.shuffle = T, condition.disjoint = T)
  
  saveRDS(each.gene.survival.res, "each.gene.survival.res.RDS")
  
  clusterEvalQ(clu, {
    each.gene.survival.res<-readRDS('each.gene.survival.res.RDS')
    NULL
  })
  file.remove('each.gene.survival.res.RDS')
  
  summary_surv_res.egfr_mutant.shuffle<-do.call(rbind, pblapply(c(1:dim(gene_dictionary)[1]), FUN=survival.marker.stat.summary, 
                                                                discovery.set=evidence.datasets, cl=clu)) %>%
    filter(n.evidence>=length(evidence.datasets)*0.8) %>%
    mutate(Previous_markers=sapply(gene,
                                   FUN=function(x){
                                     return(ifelse(!is.na(match(x, unlist(previous_markers))), "Previous marker", "Not a previous marker"))
                                   })) %>% 
    mutate(adjusted.p=p.adjust(p.value, "BH")) %>%
    mutate(significant=sapply(adjusted.p, FUN=function(x){
      return(ifelse(x<=0.05, "q \u2264 0.05", "q > 0.05"))
    }))
  
  if(length((summary_surv_res.egfr_mutant.shuffle %>% filter(adjusted.p<=0.05))$gene)>=3){
    significant.marker.nonmut.sampling<-c(significant.marker.nonmut.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(adjusted.p<=0.05))$gene))
    median.favorable.nonmut.sampling<-c(median.favorable.nonmut.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(HR<1) %>% filter(adjusted.p<=0.05))$HR))
    median.unfavorable.nonmut.sampling<-c(median.unfavorable.nonmut.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(HR>1) %>% filter(adjusted.p<=0.05))$HR))
    num.favorable.nonmut.sampling<-c(num.favorable.nonmut.sampling, dim(summary_surv_res.egfr_mutant.shuffle %>% filter(HR<1) %>% filter(adjusted.p<=0.05))[1])
    num.unfavorable.nonmut.sampling<-c(num.unfavorable.nonmut.sampling, dim(summary_surv_res.egfr_mutant.shuffle %>% filter(HR>1) %>% filter(adjusted.p<=0.05))[1])
    c.index.nonmut.sampling<-c(c.index.nonmut.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(adjusted.p<=0.05))$c.ind.mean))
    count=count+1
    cat(paste0("non-mutant sampling ", count))
  }
  
  if(count==300){
    break
  }
}

median.favorable.nonmut.sampling[which(is.na(median.favorable.nonmut.sampling))]<-1
median.unfavorable.nonmut.sampling[which(is.na(median.unfavorable.nonmut.sampling))]<-1

set.seed(123)
significant.marker.mix.sampling<-list()
median.favorable.mix.sampling<-median.unfavorable.mix.sampling<-c()
num.favorable.mix.sampling<-num.unfavorable.mix.sampling<-c()
c.index.mix.sampling<-c()
count<-1
while(TRUE){
  each.gene.survival.res<-survival.marker.stat.individual(input_query=condition, other.shuffle = T, condition.disjoint = F)
  
  saveRDS(each.gene.survival.res, "each.gene.survival.res.RDS")
  
  clusterEvalQ(clu, {
    each.gene.survival.res<-readRDS('each.gene.survival.res.RDS')
    NULL
  })
  file.remove('each.gene.survival.res.RDS')
  
  summary_surv_res.egfr_mutant.shuffle<-do.call(rbind, pblapply(c(1:dim(gene_dictionary)[1]), FUN=survival.marker.stat.summary, 
                                                                discovery.set=evidence.datasets, cl=clu)) %>%
    filter(n.evidence>=length(evidence.datasets)*0.8) %>%
    mutate(Previous_markers=sapply(gene,
                                   FUN=function(x){
                                     return(ifelse(!is.na(match(x, unlist(previous_markers))), "Previous marker", "Not a previous marker"))
                                   })) %>% 
    mutate(adjusted.p=p.adjust(p.value, "BH")) %>%
    mutate(significant=sapply(adjusted.p, FUN=function(x){
      return(ifelse(x<=0.05, "q \u2264 0.05", "q > 0.05"))
    }))
  
  if(length((summary_surv_res.egfr_mutant.shuffle %>% filter(adjusted.p<=0.05))$gene)>=3){
    significant.marker.mix.sampling<-c(significant.marker.mix.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(adjusted.p<=0.05))$gene))
    median.favorable.mix.sampling<-c(median.favorable.mix.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(HR<1) %>% filter(adjusted.p<=0.05))$HR))
    median.unfavorable.mix.sampling<-c(median.unfavorable.mix.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(HR>1) %>% filter(adjusted.p<=0.05))$HR))
    num.favorable.mix.sampling<-c(num.favorable.mix.sampling, dim(summary_surv_res.egfr_mutant.shuffle %>% filter(HR<1) %>% filter(adjusted.p<=0.05))[1])
    num.unfavorable.mix.sampling<-c(num.unfavorable.mix.sampling, dim(summary_surv_res.egfr_mutant.shuffle %>% filter(HR>1) %>% filter(adjusted.p<=0.05))[1])
    c.index.mix.sampling<-c(c.index.mix.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(adjusted.p<=0.05))$c.ind.mean))
    count=count+1
    cat(paste0("whole patient sampling ", count))
  }
  
  if(count==300){
    break
  }
  
}


median.favorable.mix.sampling[which(is.na(median.favorable.mix.sampling))]<-1
median.unfavorable.mix.sampling[which(is.na(median.unfavorable.mix.sampling))]<-1


set.seed(123)
significant.marker.mut.sampling<-list()
median.favorable.mut.sampling<-median.unfavorable.mut.sampling<-c()
num.favorable.mut.sampling<-num.unfavorable.mut.sampling<-c()
c.index.mut.sampling<-c()
count=1

while(TRUE){
  each.gene.survival.res<-survival.marker.stat.individual(input_query=condition, egfr.shuffle = T)
  
  saveRDS(each.gene.survival.res, "each.gene.survival.res.RDS")
  
  clusterEvalQ(clu, {
    each.gene.survival.res<-readRDS('each.gene.survival.res.RDS')
    NULL
  })
  file.remove('each.gene.survival.res.RDS')
  
  summary_surv_res.egfr_mutant.shuffle<-do.call(rbind, pblapply(c(1:dim(gene_dictionary)[1]), FUN=survival.marker.stat.summary, 
                                                                discovery.set=evidence.datasets, cl=clu)) %>%
    filter(n.evidence>=length(evidence.datasets)*0.8) %>%
    mutate(Previous_markers=sapply(gene,
                                   FUN=function(x){
                                     return(ifelse(!is.na(match(x, unlist(previous_markers))), "Previous marker", "Not a previous marker"))
                                   })) %>% 
    mutate(adjusted.p=p.adjust(p.value, "BH")) %>%
    mutate(significant=sapply(adjusted.p, FUN=function(x){
      return(ifelse(x<=0.05, "q \u2264 0.05", "q > 0.05"))
    }))
  
  if(length((summary_surv_res.egfr_mutant.shuffle %>% filter(adjusted.p<=0.05))$gene)>=3){
    significant.marker.mut.sampling<-c(significant.marker.mut.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(adjusted.p<=0.05))$gene))
    median.favorable.mut.sampling<-c(median.favorable.mut.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(HR<1) %>% filter(adjusted.p<=0.05))$HR))
    median.unfavorable.mut.sampling<-c(median.unfavorable.mut.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(HR>1) %>% filter(adjusted.p<=0.05))$HR))
    num.favorable.mut.sampling<-c(num.favorable.mut.sampling, dim(summary_surv_res.egfr_mutant.shuffle %>% filter(HR<1) %>% filter(adjusted.p<=0.05))[1])
    num.unfavorable.mut.sampling<-c(num.unfavorable.mut.sampling, dim(summary_surv_res.egfr_mutant.shuffle %>% filter(HR>1) %>% filter(adjusted.p<=0.05))[1])
    c.index.mut.sampling<-c(c.index.mut.sampling, list((summary_surv_res.egfr_mutant.shuffle %>% filter(adjusted.p<=0.05))$c.ind.mean))
    count=count+1
    cat(paste0("egfr-mutant sampling ", count))
  }
  
  if(count==300){
    break
  }
  
}

median.favorable.mut.sampling[which(is.na(median.favorable.mut.sampling))]<-1
median.unfavorable.mut.sampling[which(is.na(median.unfavorable.mut.sampling))]<-1

rm(list=ls()[-match(c("c.index.mix.sampling", 'c.index.nonmut.sampling', 'num.favorable.mix.sampling', 'num.favorable.nonmut.sampling',
        'num.unfavorable.mix.sampling', 'num.unfavorable.nonmut.sampling', 'median.favorable.mix.sampling', 'median.favorable.nonmut.sampling',
        'median.unfavorable.mix.sampling', 'median.unfavorable.nonmut.sampling', 'summary_surv_res.egfr_mutant',
        'significant.marker.mix.sampling', 'significant.marker.nonmut.sampling', 
        'c.index.mut.sampling', 'median.favorable.mut.sampling', 'median.unfavorable.mut.sampling', 'significant.marker.mut.sampling'), ls())])
