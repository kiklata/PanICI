library(Seurat)

CD8.Tex.harmony <- readRDS("~/PaperCD8/data/Tex/CD8.Tex.harmony.rds")

Idents(CD8.Tex.harmony) = CD8.Tex.harmony$manual.celltype.Tex

Tex.downsample = subset(CD8.Tex.harmony,downsample = 2500) 

saveRDS(Tex.downsample,file = "~/PaperCD8/data/Tex/CD8.Tex.downsample.rds")
