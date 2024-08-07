library(tibble)

gene_dictionary<-readRDS("data/gene_dictionary.RDS")
symbol.to.entrez<-function(symbol){
  gene_dictionary$ENTREZID[match(symbol, gene_dictionary$gene_symbol)]
}
entrez.to.symbol<-function(entrez){
  gene_dictionary$gene_symbol[match(entrez, gene_dictionary$ENTREZID)]
}


### MSigDB Category, Data source, Gene set name
# Hallmark
Category="H: Hallmark gene sets"
Source="MSigDB"
Annotation="Functional ontology"
file.read<-readLines("MSigDB/h.all.v2023.2.Hs.entrez.gmt")

f.msigdb<-do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1], Annotation=Annotation)
}))

# Curated gene sets
Category="C2: Curated gene sets"

## CGP: Chemical and genetic perturbations
Source="CGP: Chemical and genetic perturbations"
file.read<-readLines("MSigDB/c2.cgp.v2023.2.Hs.entrez.gmt")
Annotation="Molecular relation - Genetic"
f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1], Annotation=Annotation)
})))

## CP: Canonical Pathways
Source="CP: Canonical pathways"
file.read<-readLines("MSigDB/c2.cp.v2023.2.Hs.entrez.gmt")
Annotation="Pathway - Composite"

f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1],Annotation=Annotation)
})))

# Regulatory target gene sets
Category="C3: regulatory target gene sets"

## MIR: microRNA targets
Source="MIR: microRNA targets"
file.read<-readLines("MSigDB/c3.mir.v2023.2.Hs.entrez.gmt")
Annotation="Molecular relation - MicroRNA"

f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1], Annotation=Annotation)
})))

## TFT: transcription factor targets
Source="TFT: transcription factor targets"
file.read<-readLines("MSigDB/c3.tft.v2023.2.Hs.entrez.gmt")
Annotation="Molecular relation - TF"
f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1], Annotation=Annotation)
})))


# Computational gene sets
Category="C4: Computational gene sets"
Annotation="Molecular Relation - Genetic"
## 3CA: Curated Cancer Cell Atlas
Source="3CA: Curated Cancer Cell Atlas"
file.read<-readLines("MSigDB/c4.3ca.v2023.2.Hs.entrez.gmt")

f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1], Annotation=Annotation)
})))

## CGN: cancer gene neighborhoods
Source="CGN: cancer gene neighborhoods"
file.read<-readLines("MSigDB/c4.cgn.v2023.2.Hs.entrez.gmt")

f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1], Annotation=Annotation)
})))

## CM: cancer modules
Source="CM: cancer modules"
file.read<-readLines("MSigDB/c4.cgn.v2023.2.Hs.entrez.gmt")

f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1], Annotation=Annotation)
})))


# Ontology gene sets
Category="C5: Ontology gene sets"
## Gene Ontology
Source="GO: Gene Ontology gene sets"
file.read<-readLines("MSigDB/c5.go.v2023.2.Hs.entrez.gmt")
Annotation="Functional ontology"

f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1], Annotation=Annotation)
})))

## Human Phenotype Ontology
Source="HPO: Human Phenotype Ontology"
file.read<-readLines("MSigDB/c5.hpo.v2023.2.Hs.entrez.gmt")
Annotation="Phenotype ontology"

f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1], Annotation=Annotation)
})))

# Oncogenic signature gene sets
Category="C6: Oncogenic signature gene sets"
Source="Oncogenic signature gene sets"
file.read<-readLines("MSigDB/c6.all.v2023.2.Hs.entrez.gmt")
Annotation="Molecular relation - Genetic"

f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1], Annotation=Annotation)
})))

# Immunologic signature gene sets
Category="C7: Immunologic signature gene sets"
Annotation="Functional ontology"

Source="ImmuneSigDB"
file.read<-readLines("MSigDB/c7.immunesigdb.v2023.2.Hs.entrez.gmt")
f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(DB="MSigDB", Source=Source, Name=line[1], Annotation=Annotation)
})))


if(FALSE){
  ## VAX: vaccine response gene sets
  Source="VAX: vaccine response gene sets"
  file.read<-readLines("MSigDB/c7.vax.v2023.2.Hs.entrez.gmt")
  f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
    line=unlist(strsplit(x, "\t", fixed=T))[-2]
    tibble(DB="MSigDB", Source=Source, Name=line[1])
  })))
  ## C8: cell type signature gene sets
  Category="C8: cell type signature gene sets"
  Source="cell type signature gene sets"
  file.read<-readLines("MSigDB/c8.all.v2023.2.Hs.entrez.gmt")
  f.msigdb<-rbind(f.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
    line=unlist(strsplit(x, "\t", fixed=T))[-2]
    tibble(DB="MSigDB", Source=Source, Name=line[1])
  })))
}


f.msigdb<-as_tibble(cbind(MOD_ID=c(1:dim(f.msigdb)[1]), f.msigdb))


# Hallmark
Category="H: Hallmark gene sets"
Source="MSigDB"
file.read<-readLines("MSigDB/h.all.v2023.2.Hs.entrez.gmt")

f2g.msigdb<-do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
}))


# Curated gene sets
Category="C2: Curated gene sets"

## CGP: Chemical and genetic perturbations
Source="CGP: Chemical and genetic perturbations"
file.read<-readLines("MSigDB/c2.cgp.v2023.2.Hs.entrez.gmt")

f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
})))

## CP: Canonical Pathways
Source="CP: Canonical pathways"
file.read<-readLines("MSigDB/c2.cp.v2023.2.Hs.entrez.gmt")

f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
})))


# Regulatory target gene sets
Category="C3: regulatory target gene sets"

## MIR: microRNA targets
Source="MIR: microRNA targets"
file.read<-readLines("MSigDB/c3.mir.v2023.2.Hs.entrez.gmt")

f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
})))


## TFT: transcription factor targets
Source="TFT: transcription factor targets"
file.read<-readLines("MSigDB/c3.tft.v2023.2.Hs.entrez.gmt")

f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
})))

# Computational gene sets
Category="C4: Computational gene sets"

## 3CA: Curated Cancer Cell Atlas
Source="3CA: Curated Cancer Cell Atlas"
file.read<-readLines("MSigDB/c4.3ca.v2023.2.Hs.entrez.gmt")

f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
})))


## CGN: cancer gene neighborhoods
Source="CGN: cancer gene neighborhoods"
file.read<-readLines("MSigDB/c4.cgn.v2023.2.Hs.entrez.gmt")

f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[intersect(which(f.msigdb$Name==line[1]), which(f.msigdb$Source==Source))], Entrez=line[-1])
})))


## CM: cancer modules
Source="CM: cancer modules"
file.read<-readLines("MSigDB/c4.cgn.v2023.2.Hs.entrez.gmt")

f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[intersect(which(f.msigdb$Name==line[1]), which(f.msigdb$Source==Source))], Entrez=line[-1])
})))



# Ontology gene sets
Category="C5: Ontology gene sets"
## Gene Ontology
Source="GO: Gene Ontology gene sets"
file.read<-readLines("MSigDB/c5.go.v2023.2.Hs.entrez.gmt")
f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
})))


## Human Phenotype Ontology
Source="HPO: Human Phenotype Ontology"
file.read<-readLines("MSigDB/c5.hpo.v2023.2.Hs.entrez.gmt")

f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
})))

# Oncogenic signature gene sets
Category="C6: Oncogenic signature gene sets"
Source="Oncogenic signature gene sets"
file.read<-readLines("MSigDB/c6.all.v2023.2.Hs.entrez.gmt")
f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
})))

# Immunologic signature gene sets
Category="C7: Immunologic signature gene sets"

Source="ImmuneSigDB"
file.read<-readLines("MSigDB/c7.immunesigdb.v2023.2.Hs.entrez.gmt")
f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
  line=unlist(strsplit(x, "\t", fixed=T))[-2]
  tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
})))

if(FALSE){
  
  ## VAX: vaccine response gene sets
  Source="VAX: vaccine response gene sets"
  file.read<-readLines("MSigDB/c7.vax.v2023.2.Hs.entrez.gmt")
  f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
    line=unlist(strsplit(x, "\t", fixed=T))[-2]
    tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
  })))
  
  
  ## C8: cell type signature gene sets
  Category="C8: cell type signature gene sets"
  Source="cell type signature gene sets"
  file.read<-readLines("MSigDB/c8.all.v2023.2.Hs.entrez.gmt")
  f2g.msigdb<-rbind(f2g.msigdb, do.call(rbind, lapply(file.read, FUN=function(x){
    line=unlist(strsplit(x, "\t", fixed=T))[-2]
    tibble(MOD_ID=f.msigdb$MOD_ID[which(f.msigdb$Name==line[1])], Entrez=line[-1])
  })))
  
}

### EnrichR curation
enrichr.list<-read.csv("enrichr_list.CSV")

f.enrichr<-do.call(rbind, pblapply(c(1:dim(enrichr.list)[1]), FUN=function(i){
  url=paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", enrichr.list$Name[i])
  file.read<-readLines(con = url)
  
  res<-do.call(rbind, lapply(file.read, FUN=function(x){
    line=unlist(strsplit(x, "\t", fixed=T))[-2]
    tibble(DB="EnrichR", Source=enrichr.list$Name[i], Name=line[1], Annotation=enrichr.list$Annotation[i])
  }))
  
  return(res)
}))

f.enrichr<-as_tibble(cbind(MOD_ID=c(1:dim(f.enrichr)[1])+dim(f.msigdb)[1], f.enrichr))

f2g.enrichr<-do.call(rbind, pblapply(c(1:dim(enrichr.list)[1]), FUN=function(i){
  url=paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", enrichr.list$Name[i])
  file.read<-readLines(con = url)
  Source=enrichr.list$Name[i]
  res<-do.call(rbind, lapply(file.read, FUN=function(x){
    line=unlist(strsplit(x, "\t", fixed=T))[-2]
    tibble(MOD_ID=f.enrichr$MOD_ID[intersect(which(f.enrichr$Name==line[1]), which(f.enrichr$Source==Source))],
           Entrez=na.omit(symbol.to.entrez(line[-1])))
  }))
  return(res)
}))

f<-as_tibble(rbind(f.msigdb, f.enrichr))
f2g<-as_tibble(rbind(f2g.msigdb, f2g.enrichr))

f$mod.len<-pbsapply(f$MOD_ID, FUN=function(k){
  length(f2g[which(f2g$MOD_ID==k),]$Entrez)
}, cl=clu)

saveRDS(f, "f.RDS")
saveRDS(f2g, "f2g.RDS")

rm(Annotation);rm(Category);rm(enrichr.list);rm(f.enrichr);rm(f.msigdb)
rm(f2g.enrichr);rm(f2g.msigdb);rm(file.read);rm(Source)
