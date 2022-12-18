library(Seurat)
library(CytoTRACE)
library(dplyr)

CD8.downsample = readRDS('CD8.downsample.rds')

count =  as.matrix(CD8.downsample@assays$RNA@counts)

res <-  CytoTRACE(mat = count,ncores = 8,subsamplesize = 2000)
saveRDS(res,file = 'cytotrace.downsample.rds')

