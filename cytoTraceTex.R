library(Seurat)
library(CytoTRACE)
library(dplyr)

CD8.Tex.harmony <- readRDS("~/PaperCD8/data/Tex/CD8.Tex.downsample.rds")

count =  as.matrix(CD8.Tex.harmony@assays$RNA@counts)

res <-  CytoTRACE(mat = count,ncores = 8,subsamplesize = 2000)
saveRDS(res,file = 'cytotrace.Tex.rds')

