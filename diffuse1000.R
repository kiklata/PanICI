library(Seurat)
library(dplyr)
library(destiny)

CD8.downsample = readRDS('CD8.downsample.rds')

GeneExp.mtx <- GetAssayData(CD8.downsample, assay = "RNA", slot = "data") %>% as.matrix() # normalized data matrix

dm <- DiffusionMap(t(GeneExp.mtx),n_pcs = 50)

saveRDS(dm,file = 'CD8.downsample.DM.rds')