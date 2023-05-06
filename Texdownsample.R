library(Seurat)
CD8 <- readRDS("~/PaperCD8/data/reCD8.finished.rds")
CD8.Tex.downsample <- readRDS("~/PaperCD8/data/Tex/CD8.Tex.downsample.rds")
CD8.Tex.harmony <- readRDS("~/PaperCD8/data/Tex/CD8.Tex.harmony.rds")

tn = subset(CD8,manual.celltype.minor == 'Tn.c01.CCR7')
tn = subset(tn,downsample = 3000)
tem1 = subset(CD8,manual.celltype.minor == 'Tem.c02.GZMK' )
tem1 = subset(tem1,downsample = 2000)
tem2 = subset(CD8,manual.celltype.minor == 'Tem.c04.GZMH')
tem2 = subset(tem2,downsample = 2000)
trm = subset(CD8,manual.celltype.minor == 'Trm.c07.ZNF683')
trm = subset(trm,downsample = 1500)

t1 = merge(tn,c(tem1,tem2,trm,CD8.Tex.downsample))


CD8.Tex.harmony$label = CD8.Tex.harmony$manual.celltype.minor


saveRDS(t1,file = "~/PaperCD8/data/Tex/CD8.downsample.rds")
