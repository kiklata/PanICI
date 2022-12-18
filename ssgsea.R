library(Seurat)
library(dplyr)

obj=readRDS("~/PaperCD8/data/reCD8.finished.rds")
CD8.Tex.harmony <- readRDS("~/CD8.Tex.harmony.rds")

CD8.Tex.label = as.data.frame(CD8.Tex.harmony$label)
c1cell = rownames(filter(CD8.Tex.label,`CD8.Tex.harmony$label` == 'Tex.c1'))
c2cell = rownames(filter(CD8.Tex.label,`CD8.Tex.harmony$label` == 'Tex.c2'))
c3cell = rownames(filter(CD8.Tex.label,`CD8.Tex.harmony$label` == 'Tex.c3'))

CD8 = subset(obj,manual.celltype.major %in% c('Naive','Effector memory','MHC II','Resident memory','Exhausted'))

CD8$manual.celltype.minor = if_else(colnames(CD8) %in% c1cell,'Tex.c1',
                                    if_else(colnames(CD8) %in% c2cell,'Tex.c2',
                                            if_else(colnames(CD8) %in% c3cell,'Tex.c3',
                                                    as.character(CD8$manual.celltype.minor))))

Idents(CD8) = CD8$manual.celltype.minor

stem.feature = c('CCR7','TCF7','SELL','LEF1','CCR5')
resident.feature = c('NR4A1','NR4A3','CD69','CXCR6','ITGAE')
cyto.feature = c('IFNG','GNLY','GZMB','GZMK','GZMH','GZMA','NKG7','FGFBP2')
exhaust.feature =c('TOX2','SOX4','TIGIT','PDCD1','CTLA4','HAVCR2','LAG3','CXCL13')
costi.feature = c('ICOS','TNFSF14','TNFRSF25','TNFRSF9','CD28','TNFSF4')

sig.list=list(stem = stem.feature,resident = resident.feature,cyto = cyto.feature,exhaust = exhaust.feature,costi = costi.feature)

count = CD8@assays$RNA@counts
library(GSVA)
gsva_es <- gsva(as.matrix(count), min.sz = 3,sig.list, method=c("ssgsea"), kcdf=c("Poisson")) 
signature_exp<-as.matrix(gsva_es)
saveRDS(signature_exp,file = 'manual.ssgsea.rds')

CD8@assays$score$score = signature_exp

saveRDS(CD8,file = 'CD8.ssgsea.sig.rds')
