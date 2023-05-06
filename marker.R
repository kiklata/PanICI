library(Seurat)

library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 2400 * 1024^2)
options(future.rng.onMisuse="ignore")

obj=readRDS("~/PaperCD8/data/reCD8.finished.rds")

Idents(obj) = obj$manual.celltype.major
major.marker = FindAllMarkers(obj)
saveRDS(major.marker,file = 'major.marker.rds')

Idents(obj) = obj$manual.celltype.minor
minor.marker = FindAllMarkers(obj)
saveRDS(minor.marker,file = 'minor.marker.rds')


CD8=readRDS("~/PaperCD8/data/reCD8.finished.rds")

CD8.res = subset(CD8,treatment.efficacy == 'R')
Idents(CD8.res) = CD8.res$sample.timepoint

CD8.all.res.marker = FindAllMarkers(CD8.res)
CD8.tumor.res = subset(CD8.res,sample.Tissue == 'Tumor')
CD8.tumor.res.marker = FindAllMarkers(CD8.tumor.res)

CD8.Tex.harmony <- readRDS("~/CD8.Tex.harmony.rds")
Idents(CD8.Tex.harmony) = CD8.Tex.harmony$label
CD8.tex.res.marker = FindAllMarkers(CD8.Tex.harmony)
saveRDS(CD8.tex.res.marker,file = 'Texmarker.rds')
