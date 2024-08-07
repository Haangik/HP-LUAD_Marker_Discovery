### Functional annotation analysis
f<-readRDS("f.RDS")
f2g<-readRDS("f2g.RDS")
clusterEvalQ(clu, {
  f<-readRDS("f.RDS")
  f2g<-readRDS("f2g.RDS")
  NULL
})


total.genes<-intersect(gene_dictionary$ENTREZID, unique(f2g$Entrez))

clusterExport(clu, 'total.genes')

target.mod.id<-(f %>% filter(mod.len<=length(total.genes)*0.05))$MOD_ID

target.functions<-f %>% filter(MOD_ID %in% target.mod.id)

TableS4<-target.functions %>% dplyr::select(DB, Source, Annotation) %>% distinct()
TableS4$Num_Set<-sapply(TableS4$Source, FUN=function(x){length(which(target.functions$Source %in% x))})



enrichment.inhouse.mod_egfr.mutants<-do.call(rbind, 
                                             pblapply(target.mod.id, FUN=inhouse.annotation,
                                                      target.gene=symbol.to.entrez(na.omit(summary_surv_res.egfr_mutant %>% filter(adjusted.p<=0.05))$gene),
                                                      cl=clu)) %>% mutate(adjusted.p=p.adjust(pval))

enrichment.inhouse.mod_egfr.wt<-do.call(rbind, 
                                             pblapply(target.mod.id, FUN=inhouse.annotation,
                                                      target.gene=symbol.to.entrez(na.omit(summary_surv_res.egfr_wt %>% filter(adjusted.p<=0.05))$gene),
                                                      cl=clu)) %>% mutate(adjusted.p=p.adjust(pval))

markers.random.sampling.mut<-
  significant.marker.mut.sampling[which(sapply(significant.marker.mut.sampling, length)>=3)]

clusterExport(clu, "markers.random.sampling.mut")
enrichment.inhouse.mod_random.markers.mut<-
  lapply(c(1:length(markers.random.sampling.mut)), FUN=function(i){
    res<-do.call(rbind, 
                 pblapply(target.mod.id, FUN=inhouse.annotation,
                          target.gene=na.omit(symbol.to.entrez(markers.random.sampling.mut[[i]])),
                          cl=clu))%>% mutate(adjusted.p=p.adjust(pval)) %>% filter(adjusted.p<=0.05)
    return(res)
  })


markers.random.sampling.mix<-
  significant.marker.mix.sampling[which(sapply(significant.marker.mix.sampling, length)>=3)]

clusterExport(clu, "markers.random.sampling.mix")
enrichment.inhouse.mod_random.markers.mix<-
  lapply(c(1:length(markers.random.sampling.mix)), FUN=function(i){
    res<-do.call(rbind, 
                 pblapply(target.mod.id, FUN=inhouse.annotation,
                          target.gene=na.omit(symbol.to.entrez(markers.random.sampling.mix[[i]])),
                          cl=clu))%>% mutate(adjusted.p=p.adjust(pval)) %>% filter(adjusted.p<=0.05)
    return(res)
  })


markers.random.sampling.nonmut<-
  significant.marker.nonmut.sampling[which(sapply(significant.marker.nonmut.sampling, length)>=3)]
clusterExport(clu, "markers.random.sampling.nonmut")

enrichment.inhouse.mod_random.markers.nonmut<-
  lapply(c(1:length(markers.random.sampling.nonmut)), FUN=function(i){
    cat(paste0(i, "\n"))
    res<-do.call(rbind, 
                 pblapply(target.mod.id, FUN=inhouse.annotation,
                          target.gene=na.omit(symbol.to.entrez(markers.random.sampling.nonmut[[i]])),
                          cl=clu))%>% mutate(adjusted.p=p.adjust(pval)) %>% filter(adjusted.p<=0.05)
    return(res)
  })

egfr_mut_sig_mod<-enrichment.inhouse.mod_egfr.mutants$mod.id[which(enrichment.inhouse.mod_egfr.mutants$adjusted.p<=0.05)]

proportion=0.15

random_marker_mod_mut_freq<-(lapply(enrichment.inhouse.mod_random.markers.mut, FUN=function(x){
  return(x$mod.id)
}) %>% unlist %>% table %>% sort)

random_marker_sig_mod_mut <-
  names(random_marker_mod_mut_freq)[which(random_marker_mod_mut_freq>=length(enrichment.inhouse.mod_random.markers.mut)*proportion)] %>%
  as.numeric %>% sort



random_marker_mod_mix_freq<-(lapply(enrichment.inhouse.mod_random.markers.mix, FUN=function(x){
       return(x$mod.id)
   }) %>% unlist %>% table %>% sort)

random_marker_sig_mod_mix <- 
  names(random_marker_mod_mix_freq)[which(random_marker_mod_mix_freq>=length(enrichment.inhouse.mod_random.markers.mut)*proportion)] %>% 
  as.numeric %>% sort


random_marker_mod_nonmut_freq<-(lapply(enrichment.inhouse.mod_random.markers.nonmut, FUN=function(x){
  return(x$mod.id)
}) %>% unlist %>% table %>% sort)

random_marker_sig_mod_nonmut<-
  names(random_marker_mod_nonmut_freq)[which(random_marker_mod_nonmut_freq>=length(enrichment.inhouse.mod_random.markers.mut)*proportion)] %>% 
  as.numeric %>% sort

random.sig.mod.mut.ind<-which(f$MOD_ID %in% random_marker_sig_mod_mut)
random.sig.mod.wt.ind<-which(f$MOD_ID %in% random_marker_sig_mod_nonmut)
random.sig.mod.entire.ind<-which(f$MOD_ID %in% random_marker_sig_mod_mix)

TableS5a<-data.frame(Source=f$Source[random.sig.mod.mut.ind], 
           Database=f$DB[random.sig.mod.mut.ind],
           Name=f$Name[random.sig.mod.mut.ind],
           Category=f$Annotation[random.sig.mod.mut.ind],
           n.genes=f$mod.len[random.sig.mod.mut.ind],
           Count=as.numeric(random_marker_mod_mut_freq)[match(random_marker_sig_mod_mut, 
                                                  names(random_marker_mod_mut_freq))]) %>%
  arrange(-Count)

TableS5b<-data.frame(Source=f$Source[random.sig.mod.wt.ind], 
                     Database=f$DB[random.sig.mod.wt.ind],
                     Name=f$Name[random.sig.mod.wt.ind],
                     Category=f$Annotation[random.sig.mod.wt.ind],
                     n.genes=f$mod.len[random.sig.mod.wt.ind],
                     Count=as.numeric(random_marker_mod_nonmut_freq)[match(random_marker_sig_mod_nonmut, 
                                                                        names(random_marker_mod_nonmut_freq))]) %>%
  arrange(-Count)

TableS5c<-data.frame(Source=f$Source[random.sig.mod.entire.ind], 
                     Database=f$DB[random.sig.mod.entire.ind],
                     Name=f$Name[random.sig.mod.entire.ind],
                     Category=f$Annotation[random.sig.mod.entire.ind],
                     n.genes=f$mod.len[random.sig.mod.entire.ind],
                     Count=as.numeric(random_marker_mod_mix_freq)[match(random_marker_sig_mod_mix, 
                                                                           names(random_marker_mod_mix_freq))]) %>%
  arrange(-Count)


write.csv(TableS5a, "TableS5a.csv", row.names = F)
write.csv(TableS5b, "TableS5b.csv", row.names = F)
write.csv(TableS5c, "TableS5c.csv", row.names = F)


target.module.id<-unique(c(random_marker_sig_mod_mut, 
                           random_marker_sig_mod_mix,
                           random_marker_sig_mod_nonmut)) %>% sort

module_gene_mat<-do.call(rbind, pblapply(target.module.id, FUN=function(i){
  filtered.mod<-f2g[which(f2g$MOD_ID==i), ]
  vec<-rep(0, length(total.genes))
  names(vec)<-total.genes
  modgene<-intersect(filtered.mod$Entrez, total.genes)
  vec[match(modgene, names(vec))]<-1
  return(vec)
}, cl=clu))

remove.ind<-which(colSums(module_gene_mat)==0)
module_gene_mat<-module_gene_mat[,-remove.ind]
total.genes.sig.mod<-total.genes[-remove.ind]
  

data.umap<-umap2(module_gene_mat, metric="cosine", n_threads=25, seed=123, n_sgd_threads = 1)
df.module.cluster<-data.frame(mod_id=target.module.id, data.umap,
                              label=factor(pbsapply(target.module.id, FUN=function(i){
                                if(i %in% egfr_mut_sig_mod){
                                  return("EGFR-mutated")
                                }else if(i %in% random_marker_sig_mod_nonmut){
                                  res<-"Non-EGFR-mutated"
                                }else{
                                  res<-"Entire patients"
                                }
                                if(i %in% intersect(random_marker_sig_mod_nonmut, random_marker_sig_mod_mix)){
                                  res<-sample(c("Non-EGFR-mutated", "Entire patients"), 1)
                                }
                                
                                return(res)
                              }), levels = c("Entire patients",
                                             "Non-EGFR-mutated", 
                                             "EGFR-mutated"
                              )),
                              mod_len=sapply(target.module.id, FUN=function(k){
                                f$mod.len[which(f$MOD_ID==k)]
                              })
)

params<-do.call(rbind, lapply(seq(0.1, 0.5, 0.05), FUN=function(eps){
  data.frame(eps=eps, minPts=seq(30, 150, 5))
}))
clusterExport(clu, c("df.module.cluster","params"))


dbcv.tuning<-do.call(rbind, pblapply(c(1:dim(params)[1]), FUN=function(k){
  randn<-sample(c(1:1000000),1)
  eps=params$eps[k]
  minPts=params$minPts[k]
  x<-as.matrix(df.module.cluster[,c(2,3)])
  db<-dbscan(x, eps=eps, minPts=minPts)
  df.module.cluster$cluster<-db$cluster
  write.table(df.module.cluster %>% dplyr::select(X1, X2),
              file = paste0("module_cluster_embed_", randn, ".txt"), quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(df.module.cluster %>% dplyr::select(cluster),
              file = paste0("module_cluster_label_", randn,".txt"), quote = F, sep = "\t", row.names = F, col.names = F)
  dbcv.score<-dbcv_execute(randn)
  file.remove(paste0("module_cluster_embed_", randn, ".txt"))
  file.remove(paste0("module_cluster_label_", randn,".txt"))
  
  py_gc <- import("gc")
  py_gc$collect()
  return(data.frame(eps, minPts, dbcv.score))
}, cl=clu)) %>% arrange(-dbcv.score)



eps=dbcv.tuning$eps[1]
minPts=dbcv.tuning$minPts[1]
db<-dbscan(as.matrix(df.module.cluster[,c(2,3)]), eps=eps, minPts=minPts)
df.module.cluster$cluster<-sapply(db$cluster, FUN=function(x){
  if(x!=0){
    factor(paste0("Cluster ", x, collapse = ""),
           levels=c(sapply(setdiff(names(table(db$cluster)), "0"), FUN=function(num){paste0("Cluster ", num)}),
                    "Outliers"))  
  }else{
    factor("Outliers", 
           levels=c(sapply(setdiff(names(table(db$cluster)), "0"), FUN=function(num){paste0("Cluster ", num)}),
                    "Outliers"))
  }
})


Figure3a<-ggplot(df.module.cluster %>% arrange(label))+
  theme_classic()+
  geom_point(aes(x=X1, y=X2, color=label), size=1)+
  xlab("")+
  ylab("")+
  ggtitle("Clustering of significant functional modules based on gene membership similarity")+
  scale_color_manual(breaks=c("EGFR-mutated", "Non-EGFR-mutated", "Entire patients"), values = c( 'red',"#619cff", "black"))+
  labs(color = "EGFR mutation status")+
  geom_mark_hull(data=df.module.cluster %>% arrange(label) %>% filter(cluster!="Outliers"),
                 aes(x=X1, y=X2, group=cluster, label=cluster), concavity=6, show.legend = FALSE,
                 label.fontsize=10,
                 con.size=0.2,con.type="straight", con.border = 'none', 
                 expand=unit(1.3, "mm"),  linetype="dotted", size=0.25)+
  theme(text=element_text(size=11), legend.key.size = unit(0.1, 'cm'))


tiff(file="Figure3a.tiff", width = 20, height = 17, units = "cm", res=300)
Figure3a
dev.off()

egfr_mutant_markers_mod_cluster_perc<-
  df.module.cluster$cluster[na.omit(match(egfr_mut_sig_mod, df.module.cluster$mod_id))] %>% table

egfr_mutant_markers_mod_cluster_perc<-
  data.frame(Percentage=100*as.numeric(egfr_mutant_markers_mod_cluster_perc/sum(egfr_mutant_markers_mod_cluster_perc)),
             Cluster=names(egfr_mutant_markers_mod_cluster_perc),
             Label="EGFR-mutated")

random_patient_markers_mut_mod_cluster_perc<-do.call(rbind, lapply(enrichment.inhouse.mod_random.markers.mut, FUN=function(x){
  tb<-table(df.module.cluster$cluster[na.omit(match(x$mod.id, df.module.cluster$mod_id))])
  if(sum(tb)!=0)
    tb/sum(tb)
  else
    return()
}))

random_patient_markers_mut_mod_cluster_perc<-colMeans(random_patient_markers_mut_mod_cluster_perc)
random_patient_markers_mut_mod_cluster_perc<-data.frame(Percentage=100*random_patient_markers_mut_mod_cluster_perc, 
                                                        Cluster=names(random_patient_markers_mut_mod_cluster_perc),
                                                        Label="EGFR-mutated")

random_patient_markers_mix_mod_cluster_perc<-do.call(rbind, lapply(enrichment.inhouse.mod_random.markers.mix, FUN=function(x){
  tb<-table(df.module.cluster$cluster[na.omit(match(x$mod.id, df.module.cluster$mod_id))])
  if(sum(tb)!=0)
    tb/sum(tb)
  else
    return()
}))

random_patient_markers_mix_mod_cluster_perc<-colMeans(random_patient_markers_mix_mod_cluster_perc)
random_patient_markers_mix_mod_cluster_perc<-data.frame(Percentage=100*random_patient_markers_mix_mod_cluster_perc, 
                                          Cluster=names(random_patient_markers_mix_mod_cluster_perc),
                                          Label="Entire patients")

random_patient_markers_nonmut_mod_cluster_perc<-do.call(rbind, lapply(enrichment.inhouse.mod_random.markers.nonmut, FUN=function(x){
  tb<-table(df.module.cluster$cluster[na.omit(match(x$mod.id, df.module.cluster$mod_id))])
  if(sum(tb)!=0)
    tb/sum(tb)
  else
    return()
  
}))

random_patient_markers_nonmut_mod_cluster_perc<-colMeans(random_patient_markers_nonmut_mod_cluster_perc)
random_patient_markers_nonmut_mod_cluster_perc<-data.frame(Percentage=100*random_patient_markers_nonmut_mod_cluster_perc, 
                                                        Cluster=names(random_patient_markers_nonmut_mod_cluster_perc),
                                                        Label="Non-EGFR-mutated")

df.percentage<-rbind(#egfr_mut_sig_mod,
  random_patient_markers_mut_mod_cluster_perc,
                     random_patient_markers_mix_mod_cluster_perc,
                     random_patient_markers_nonmut_mod_cluster_perc)

other_clusters<-df.percentage %>% filter(!Cluster %in% c("Cluster 1", "Cluster 3", "Outliers"))


df.percentage<-rbind(df.percentage %>% filter(Cluster %in% c("Cluster 1", "Cluster 3", "Outliers")), 
                     data.frame(Percentage=sapply(unique(other_clusters$Label), FUN=function(x){
                       (other_clusters %>% filter(Label==x))$Percentage %>% sum}),
                       Cluster="Other Clusters",
                       Label=unique(other_clusters$Label)))

df.percentage$Label<-factor(df.percentage$Label, levels=c("EGFR-mutated", "Non-EGFR-mutated", "Entire patients"))


Figure3b<-ggplot(df.percentage, aes(x=Cluster, y=Percentage, fill=Label))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
  labs(fill="EGFR mutation status")+
  xlab("")+ylab("Percentage (%)")+
  scale_y_continuous(limits = c(0, 100))+
  scale_fill_manual(values=c("red", "#619cff", "black" ))+
  theme(text=element_text(size=13))

tiff(file="Figure3b.tiff", width = 18, height = 12, units = "cm", res=300)
Figure3b
dev.off()

modules_per_clusters<-lapply(levels(df.module.cluster$cluster), FUN=function(c){
  data.frame(Cluster=c, f[match((df.module.cluster %>% filter(cluster==c))$mod_id, f$MOD_ID) %>% sort, ])
})

write.table(do.call(rbind, modules_per_clusters), file = "Modules_per_cluster.txt", row.names = F, sep="\t")

clusterExport(clu, c("egfr.significant.gene.info", "egfr.wt.significant.gene.info", "whole.patients.significant.gene.info",
                     "symbol.to.entrez", "inhouse.annotation"))
cluster_module_OR<-lapply(c(1:3), FUN=function(i){
  if(i!=2){
    Cluster<-paste0(c("Cluster ", i), collapse = "")
    module<-(df.module.cluster %>% filter(cluster %in% Cluster))$mod_id
  }else{
    module<-(df.module.cluster %>% filter(!cluster %in% c("Cluster 1", "Cluster 3", "Outliers")))$mod_id
  }
  module.OR<-do.call(rbind, pblapply(module, FUN=function(mod_id){
    mod_gene<-(f2g %>% filter(MOD_ID==mod_id))$Entrez
    egfr_mut_marker_OR=inhouse.annotation(n=mod_id, 
                                          target.gene =egfr.significant.gene.info$gene %>%
                                            symbol.to.entrez%>% na.omit)$OR
    egfr_wt_marker_OR=inhouse.annotation(n=mod_id, 
                                         target.gene =egfr.wt.significant.gene.info$gene %>%
                                           symbol.to.entrez%>% na.omit)$OR
    whole_marker_OR=inhouse.annotation(n=mod_id, 
                                       target.gene =whole.patients.significant.gene.info$gene %>%
                                         symbol.to.entrez%>% na.omit)$OR
    return(data.frame(egfr_mut=egfr_mut_marker_OR,
                      egfr_wt=egfr_wt_marker_OR,
                      whole=whole_marker_OR))
  }, cl=clu))
  return(module.OR)
})

df.figure3c<-data.frame(Group=c(rep("EGFR-mutated", 3), rep("Non-EGFR-mutated",3)),
                        Cluster=rep(c("Cluster 1", "Cluster 3", "Other Clusters"),2),
           OR=c(sapply(cluster_module_OR, FUN=function(x){return(median(x$egfr_mut))}),
                sapply(cluster_module_OR, FUN=function(x){return(median(x$egfr_wt))})))

Figure3c<-ggplot(df.figure3c, aes(x=Cluster, y=OR, fill=Group))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
  labs(fill="EGFR mutation status")+
  geom_hline(yintercept = 1, linetype="dotted")+
  xlab("")+ylab("Odds ratio")+
  scale_fill_manual(values=c("red", "#619cff"))+
  theme(text=element_text(size=13))

tiff(file="Figure3c.tiff", width = 15, height = 12, units = "cm", res=300)
Figure3c
dev.off()

