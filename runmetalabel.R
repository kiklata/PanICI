library(Seurat)
library(dplyr)

ALL.info.kegg.remove.blood <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/obj/ALL.info.kegg.remove.blood.rds")
before = subset(ALL.info.kegg.remove.blood,sample.timepoint == 'Before')

ALL.info.hgnc <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/obj/ALL.info.hgnc.rds")
BC_Bassez.hgnc <- readRDS("~/BC_Bassez.hgnc.rds")
allmarker <- readRDS("~/allmarker.rds")

BC_Zhang = subset(before,Study == 'BC_Zhang')
BC_Zhang = ALL.info.hgnc[,colnames(BC_Zhang)]

BCC = subset(before,Study == 'BCC_Yost')
BCC = ALL.info.hgnc[,colnames(BCC)]

SCC = subset(before,Study == 'SCC_Yost')
SCC = ALL.info.hgnc[,colnames(SCC)]

genebl = readRDS('genebl.rds')

final.gene = setdiff(rownames(BC_Zhang),y = genebl)
BC_Zhang = BC_Zhang[final.gene,]

final.gene = setdiff(rownames(BCC),y = genebl)
BCC = BCC[final.gene,]

final.gene = setdiff(rownames(SCC),y = genebl)
SCC = SCC[final.gene,]

run.metalabel = function(obj,ref){
  
  library(Seurat)
  library(dplyr)
  ref = ref %>% NormalizeData() %>% FindVariableFeatures()
  query = obj %>% NormalizeData() %>% FindVariableFeatures()
  anchors = FindTransferAnchors(reference = ref, query = query, dims = 1:30) 
  predictions <- TransferData(anchorset = anchors, refdata = ref$MetaState, dims = 1:30) 
  
  return(predictions)
}

BC_Zhang.trans = run.metalabel(obj = BC_Zhang,ref = BC_Bassez.hgnc)
BCC.trans = run.metalabel(obj = BCC,ref = BC_Bassez.hgnc)
SCC.trans = run.metalabel(obj = SCC,ref = BC_Bassez.hgnc)
save(BC_Zhang.trans,BCC.trans,SCC.trans,file = 'trans.Rdata')