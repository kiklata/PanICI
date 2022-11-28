library(Seurat)
library(dplyr)
library(harmony)

library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 2400 * 1024^2)
options(future.rng.onMisuse="ignore")

reCD8 <- readRDS("~/PaperCD8/data/reCD8.finished.rds")

Tex = subset(reCD8,manual.celltype.major == 'Exhausted')

Tex = FindNeighbors(Tex, reduction = "harmony", dims = 1:30)

Tex = FindClusters(Tex,resolution = seq(0.2,0.6,0.1))

saveRDS(Tex,file = 'CD8.Tex.harmony.rds')