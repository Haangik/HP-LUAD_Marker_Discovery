

## Univariate survival analysis across discovery datasets
### whole patients (no specific condition)
each.gene.survival.res<-survival.marker.stat.individual(input_query="")

clusterExport(clu, 'each.gene.survival.res')
clusterExport(clu, 'egfr.datasets')

### removing independent validation dataset (TCGA-LUAD) from marker discovery
evidence.datasets<-setdiff(egfr.datasets, "TCGA-LUAD")
summary_surv_res.whole_patients<-do.call(rbind, pblapply(c(1:dim(gene_dictionary)[1]), FUN=survival.marker.stat.summary, 
                                          discovery.set=evidence.datasets, cl=clu))  %>%
  filter(n.evidence>=length(evidence.datasets)*0.8) %>%
  mutate(Previous_markers=sapply(gene,
                                 FUN=function(x){
                                   return(ifelse(!is.na(match(x, unlist(previous_markers))), "Previous marker", "Not a previous marker"))
                                 })) %>% 
  mutate(adjusted.p=p.adjust(p.value, "BH"))

summary_surv_res.whole_patients <- 
  summary_surv_res.whole_patients %>%
  mutate(significant=sapply(adjusted.p, FUN=function(x){
    return(ifelse(x<=0.05, 
                  paste0("q \u2264 0.05 (", length(which(summary_surv_res.whole_patients$adjusted.p<=0.05)), " genes)"),
                  paste0("q > 0.05 (", length(which(summary_surv_res.whole_patients$adjusted.p>0.05)), " genes)")
    ))
  }))

whole.patients.significant.gene.info<-summary_surv_res.whole_patients %>% filter(adjusted.p<=0.05)


### Non-EGFR-mutated patients
condition<-"egfr_mut:EGFR-WT"
each.gene.survival.res<-survival.marker.stat.individual(input_query=condition)
clusterExport(clu, 'each.gene.survival.res')

summary_surv_res.egfr_wt<-do.call(rbind, pblapply(c(1:dim(gene_dictionary)[1]), FUN=survival.marker.stat.summary, 
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

summary_surv_res.egfr_wt <- 
  summary_surv_res.egfr_wt %>%
  mutate(significant=sapply(adjusted.p, FUN=function(x){
    return(ifelse(x<=0.05, 
                  paste0("q \u2264 0.05 (", length(which(summary_surv_res.egfr_wt$adjusted.p<=0.05)), " genes)"),
                  paste0("q > 0.05 (", length(which(summary_surv_res.egfr_wt$adjusted.p>0.05)), " genes)")
    ))
  }))
egfr.wt.significant.gene.info<- summary_surv_res.egfr_wt %>% filter(adjusted.p<0.05)

### EGFR-mutated patients
condition<-"egfr_mut:EGFR-Mut"
each.gene.survival.res<-survival.marker.stat.individual(input_query=condition)
clusterExport(clu, 'each.gene.survival.res')

summary_surv_res.egfr_mutant<-do.call(rbind, pblapply(c(1:dim(gene_dictionary)[1]), FUN=survival.marker.stat.summary, 
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

summary_surv_res.egfr_mutant$Previous_markers<-factor(summary_surv_res.egfr_mutant$Previous_markers, 
                                                      levels = unique(summary_surv_res.egfr_mutant$Previous_markers)[c(2,1)])
summary_surv_res.egfr_mutant <- 
  summary_surv_res.egfr_mutant %>%
  mutate(significant=sapply(adjusted.p, FUN=function(x){
    return(ifelse(x<=0.05, 
                  paste0("q \u2264 0.05 (", length(which(summary_surv_res.egfr_mutant$adjusted.p<=0.05)), " genes)"),
                  paste0("q > 0.05 (", length(which(summary_surv_res.egfr_mutant$adjusted.p>0.05)), " genes)")
    ))
  }))

egfr.significant.gene.info<- summary_surv_res.egfr_mutant %>% filter(adjusted.p<0.05)

# Figure 2. Combined univariate survival statistics of single gene expression across discovery datasets 
# Figure 2a
Figure2a3<-ggplot(summary_surv_res.whole_patients, aes(x=log(HR), y=-log(adjusted.p),  color=significant))+
  geom_point(size=0.1)+
  scale_color_manual(values = c("black", "lightgray"),
                     breaks = c((summary_surv_res.whole_patients$significant %>% unique)[2],
                                (summary_surv_res.whole_patients$significant %>% unique)[1]))+
  labs(color=NULL)+
  geom_hline(yintercept = -log(0.05), linetype="dashed")+
  theme_classic()+
  ggtitle("Entire patients")+
  xlim(-1.85, 1.85)+
  ylim(0,21)+
  theme(text=element_text(size=6), legend.key.size = unit(0.1, 'cm'))

Figure2a1<-ggplot(summary_surv_res.egfr_mutant, aes(x=log(HR), y=-log(adjusted.p),  color=significant))+
  geom_point(size=0.1)+
  scale_color_manual(values = c("black", "lightgray"),
                     breaks = c((summary_surv_res.egfr_mutant$significant %>% unique)[2],
                                (summary_surv_res.egfr_mutant$significant %>% unique)[1]))+
  labs(color=NULL)+
  geom_hline(yintercept = -log(0.05), linetype="dashed")+
  theme_classic()+
  ggtitle("EGFR-mutated patients")+
  xlim(-1.85, 1.85)+
  ylim(0,13)+
  theme(text=element_text(size=6), legend.key.size = unit(0.1, 'cm'), legend.key.width = unit(0.1, 'cm'))


Figure2a2<-ggplot(summary_surv_res.egfr_wt, aes(x=log(HR), y=-log(adjusted.p),  color=significant))+
  geom_point(size=0.1)+
  scale_color_manual(values = c("black", "lightgray"),
                     breaks = c((summary_surv_res.egfr_wt$significant %>% unique)[2],
                                (summary_surv_res.egfr_wt$significant %>% unique)[1]))+
  labs(color=NULL)+
  geom_hline(yintercept = -log(0.05), linetype="dashed")+
  theme_classic()+
  ggtitle("Non-EGFR-mutated patients")+
  xlim(-1.85, 1.85)+
  ylim(0,13)+
  theme(text=element_text(size=6), legend.key.size = unit(0.1, 'cm'), legend.key.width = unit(0.1, 'cm'))

tiff(file="Figure2a.tiff", width = 19, height = 4.725, units = "cm", res=300)
grid.arrange(Figure2a1, Figure2a2, Figure2a3, ncol=3)
dev.off()

## Supplementary Table 2
TableS2a<-summary_surv_res.egfr_mutant %>% 
  filter(adjusted.p<=0.05) %>% dplyr::select(gene, adjusted.p, HR, lower, upper, c.ind.min, pval.Q) %>% arrange(-c.ind.min)

TableS2b<-summary_surv_res.egfr_wt  %>% 
  filter(adjusted.p<=0.05) %>% dplyr::select(gene, adjusted.p, HR, lower, upper, c.ind.min, pval.Q) %>% arrange(-c.ind.min)

TableS2c<-summary_surv_res.whole_patients  %>% 
  filter(adjusted.p<=0.05) %>% dplyr::select(gene, adjusted.p, HR, lower, upper, c.ind.min, pval.Q) %>% arrange(-c.ind.min)

write.csv(TableS2a, "TableS2a.csv", row.names = F)
write.csv(TableS2b, "TableS2b.csv", row.names = F)
write.csv(TableS2c, "TableS2c.csv", row.names = F)

## Figure 2b
df.Figure2b2<-data.frame(HR=c((summary_surv_res.egfr_mutant %>% filter(adjusted.p<0.05) %>% filter(HR>1))$HR,
                    (summary_surv_res.egfr_wt %>% filter(adjusted.p<0.05) %>% filter(HR>1))$HR,
                    (summary_surv_res.whole_patients %>% filter(adjusted.p<0.05) %>% filter(HR>1))$HR),
               label=factor(c(rep("EGFR-Mutated", length((summary_surv_res.egfr_mutant %>% filter(adjusted.p<0.05) %>% filter(HR>1))$HR)),
                       rep("Non-EGFR-Mutated", length((summary_surv_res.egfr_wt %>% filter(adjusted.p<0.05) %>% filter(HR>1))$HR)),
                       rep("Entire patients", length((summary_surv_res.whole_patients %>% filter(adjusted.p<0.05) %>% filter(HR>1))$HR))),
                       levels=c("EGFR-Mutated", "Non-EGFR-Mutated", "Entire patients") %>% rev))

df.Figure2b1<-data.frame(HR=c((summary_surv_res.egfr_mutant %>% filter(adjusted.p<0.05) %>% filter(HR<1))$HR,
                    (summary_surv_res.egfr_wt %>% filter(adjusted.p<0.05) %>% filter(HR<1))$HR,
                    (summary_surv_res.whole_patients %>% filter(adjusted.p<0.05) %>% filter(HR<1))$HR),
               label=factor(c(rep("EGFR-Mutated", length((summary_surv_res.egfr_mutant %>% filter(adjusted.p<0.05) %>% filter(HR<1))$HR)),
                              rep("Non-EGFR-Mutated", length((summary_surv_res.egfr_wt %>% filter(adjusted.p<0.05) %>% filter(HR<1))$HR)),
                              rep("Entire patients", length((summary_surv_res.whole_patients %>% filter(adjusted.p<0.05) %>% filter(HR<1))$HR))),
                            levels=c("EGFR-Mutated", "Non-EGFR-Mutated", "Entire patients") %>% rev))


Figure2b<-ggplot(df.Figure2b1, aes(x=HR, y=label))+
  geom_boxplot(notch = T)+
  theme_classic()+
  stat_compare_means(df.Figure2b1, comparisons = list(c("EGFR-Mutated", "Non-EGFR-Mutated")),
                     label.y=-0.75, tip.length = -0.03, vjust=4, size=3, method.args = list(alternative = "less"))+
  stat_compare_means(df.Figure2b1, comparisons = list(c("EGFR-Mutated", "Entire patients")),
                     label.y=-0.95, tip.length = -0.03, vjust=4, size=3, method.args = list(alternative = "less"))+
  ggtitle("Hazard ratio distribution of the significant genes")+
  xlab("HR")+
  ylab("")+
  geom_boxplot(df.Figure2b2, mapping=aes(x=HR, y=label), notch = T)+
  stat_compare_means(data=df.Figure2b2, 
                     comparisons = list(c("EGFR-Mutated", "Non-EGFR-Mutated"), c("EGFR-Mutated", "Entire patients")),
                     size=3, method.args = list(alternative = "greater"))+
  scale_x_log10()+
  geom_vline(xintercept=1, linetype="dashed")+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=9),
        plot.title = element_text(size = 8.2))

tiff(file="Figure2b.tiff", width = 10.2, height = 8.1, units = "cm", res=300)
Figure2b
dev.off()

## Figure 2c
df.Figure2c<-data.frame(c.index=c((summary_surv_res.egfr_mutant %>% filter(adjusted.p<0.05))$c.ind.mean,
                         (summary_surv_res.egfr_wt %>% filter(adjusted.p<0.05))$c.ind.mean,
                         (summary_surv_res.whole_patients %>% filter(adjusted.p<0.05))$c.ind.mean),
               label=factor(c(rep("EGFR-Mutated", length((summary_surv_res.egfr_mutant %>% filter(adjusted.p<0.05))$c.ind.mean)),
                       rep("Non-EGFR-Mutated", length((summary_surv_res.egfr_wt %>% filter(adjusted.p<0.05))$c.ind.mean)),
                       rep("Entire patients", length((summary_surv_res.whole_patients %>% filter(adjusted.p<0.05))$c.ind.mean))),
                       levels=c("EGFR-Mutated", "Non-EGFR-Mutated", "Entire patients")))

Figure2c<-ggplot(df.Figure2c, aes(x=label, y=c.index))+
  geom_boxplot(notch = T)+
  theme_classic()+
  stat_compare_means(comparisons = list(c("EGFR-Mutated", "Non-EGFR-Mutated"), c("EGFR-Mutated", "Entire patients")), size=3)+
  ggtitle("Concordance indices of significant genes")+
  xlab("")+
  ylab("Average c-index")+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=9),
        plot.title = element_text(size = 8.2), 
        axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))

tiff(file="Figure2c.tiff", width = 6.8, height = 8.1, units = "cm", res=300)
Figure2c
dev.off()

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
