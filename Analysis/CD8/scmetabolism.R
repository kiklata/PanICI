
library(Seurat)
library(scMetabolism)
library(ggplot2)
library(rsvd)

ALL.info.hgnc <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/ALL.info.hgnc.rds")

ALL.info.hgnc<-sc.metabolism.Seurat(obj = ALL.info.hgnc, method = "ssGSEA", 
                                        imputation = F, ncores = 8, metabolism.type = "KEGG")
matrix <- ALL.info.hgnc@assays$METABOLISM$score

saveRDS(matrix,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/score.hngc.rds')