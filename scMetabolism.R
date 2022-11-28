library(Seurat)
library(scMetabolism)
library(ggplot2)
library(rsvd)

seu = readRDS('obj/all.T.filter.mt.geneblacklist.hsp.minC3F100.kegg.rds')
seu<-sc.metabolism.Seurat(obj = seu, method = "ssGSEA", 
                                      imputation = F, ncores = 1, metabolism.type = "KEGG")
matrix <- seu@assays$METABOLISM$score
saveRDS(matrix,file = 'score.all.T.rds')

