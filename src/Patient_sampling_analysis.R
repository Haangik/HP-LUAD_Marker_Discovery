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

## Figure S1
df.FigureS1a_favorable<-
  data.frame(HR=c(unlist(median.favorable.mix.sampling), 
                  unlist(median.favorable.nonmut.sampling),
                  unlist(median.favorable.mut.sampling)),
             group=c(rep("Entire patients",length(unlist(median.favorable.mix.sampling))),
                     rep("Non-EGFR-Mutated",length(unlist(median.favorable.nonmut.sampling))),
                     rep("EGFR-Mutated", length(unlist(median.favorable.mut.sampling)))))

df.FigureS1a_unfavorable<-
  data.frame(HR=c(unlist(median.unfavorable.mix.sampling), 
                  unlist(median.unfavorable.nonmut.sampling),
                  unlist(median.unfavorable.mut.sampling)),
             group=c(rep("Entire patients",length(unlist(median.unfavorable.mix.sampling))),
                     rep("Non-EGFR-Mutated",length(unlist(median.unfavorable.nonmut.sampling))),
                     rep("EGFR-Mutated", length(unlist(median.unfavorable.mut.sampling)))))

FigureS1a<-ggplot(df.FigureS1a_favorable %>% filter(HR>0.2),
       aes(x=HR, 
           y=factor(group, levels = c("Entire patients", "Non-EGFR-Mutated","EGFR-Mutated"))))+
  geom_boxplot(notch = T)+
  theme_classic()+
  stat_compare_means(df.FigureS1a_favorable, comparisons = list(c("EGFR-Mutated", "Non-EGFR-Mutated")),
                     label.y=-0.83,
                     tip.length = -0.03, vjust=4, method.args = list(alternative = "less"))+
  stat_compare_means(df.FigureS1a_favorable, comparisons = list(c("EGFR-Mutated", "Entire patients")),
                     label.y=-0.91, tip.length = -0.03, vjust=4, method.args = list(alternative = "less"))+
  ggtitle("Hazard ratios of the significant genes (Resampling)")+
  xlab("HR")+
  ylab("")+
  geom_boxplot(df.FigureS1a_unfavorable %>% filter(HR<2.5) , 
               mapping=aes(x=HR, 
                           y=factor(group, levels = c("Entire patients", "Non-EGFR-Mutated","EGFR-Mutated"))), 
               notch = T)+
  stat_compare_means(data=df.FigureS1a_unfavorable %>% filter(HR<2.5),
                     comparisons = list(c("EGFR-Mutated", "Non-EGFR-Mutated"), c("EGFR-Mutated", "Entire patients")),
                     method.args = list(alternative = "greater"))+
  scale_x_log10()+
  geom_vline(xintercept=1, linetype="dashed")+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12),
        plot.title = element_text(size = 8.2))

tiff(file="FigureS1a.tiff", width = 20, height = 10, units = "cm", res=300)
FigureS1a
dev.off()


df.Figures1b<-
  data.frame(C.index=c(unlist(c.index.mix.sampling),
                       unlist(c.index.nonmut.sampling),
                       unlist(c.index.mut.sampling)),
             group=c(rep("Entire patients",length(unlist(c.index.mix.sampling))),
                     rep("Non-EGFR-Mutated",length(unlist(c.index.nonmut.sampling))),
                     rep("EGFR-Mutated", length(unlist(c.index.mut.sampling)))))


FigureS1b<-ggplot(df.Figures1b, aes(x=group, y=C.index))+
  geom_boxplot(notch = T)+
  theme_classic()+
  stat_compare_means(comparisons = list(c("EGFR-Mutated", "Non-EGFR-Mutated"), c("EGFR-Mutated", "Entire patients")), size=6,
                     method.args = list(alternative = "greater"))+
  ggtitle("C-indices of significant genes (Resampling)")+
  xlab("")+
  ylab("Average c-index")+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12),
        plot.title = element_text(size = 15))

tiff(file="FigureS1b.tiff", width = 13, height = 14, units = "cm", res=300)
FigureS1b
dev.off()
