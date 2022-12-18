library(Seurat)
library(dplyr)
library(harmony)

library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 2400 * 1024^2)
options(future.rng.onMisuse="ignore")

library(scales)
library(RColorBrewer)

# CD8 Tex reduction -------------------------------------------------------

reCD8 <- readRDS("~/PaperCD8/data/reCD8.finished.rds")

CD8.Tex.harmony = subset(reCD8,manual.celltype.major == 'Exhausted')

CD8.Tex.harmony = CD8.Tex.harmony %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData(vars.to.regress = c('percent.mt','S.Score','G2M.Score')) %>% RunPCA() 

CD8.Tex.harmony =  RunHarmony(CD8.Tex.harmony,group.by.vars = 'sample.ID',assay.use='RNA', plot_convergence = F,theta = 4, 
                       kmeans_init_nstart=20, kmeans_init_iter_max=5000)

CD8.Tex.harmony = RunUMAP(CD8.Tex.harmony,reduction = 'harmony',dims = 1:20,reduction.name = 'umap_harmony')

CD8.Tex.harmony = FindNeighbors(CD8.Tex.harmony, reduction = "harmony", dims = 1:8)

CD8.Tex.harmony = FindClusters(CD8.Tex.harmony,resolution = 0.2)

Idents(CD8.Tex.harmony) = CD8.Tex.harmony$RNA_snn_res.0.2

CD8.Tex.harmony$RNA_snn_res.2 = NULL
CD8.Tex.harmony@meta.data[,c(25:33)] = NULL
CD8.Tex.harmony$seurat_clusters = CD8.Tex.harmony$RNA_snn_res.0.2

saveRDS(CD8.Tex.harmony,file = 'CD8.Tex.harmony.rds')

CD8.Tex.harmony$label = if_else(CD8.Tex.harmony$RNA_snn_res.0.2 == '0','Tex.c1',
                               if_else(CD8.Tex.harmony$RNA_snn_res.0.2 == '1','Tex.c2',
                                       if_else(CD8.Tex.harmony$RNA_snn_res.0.2 == '2','Tex.c3','Tex.c4')))

CD8.Tex.harmony$manual.celltype.Tex = if_else(CD8.Tex.harmony$label == 'Tex.c1','Tex.c01.CCL4',
                                              if_else(CD8.Tex.harmony$label == 'Tex.c2','Tex.c02.GZMH',
                                                      if_else(CD8.Tex.harmony$label == 'Tex.c3','Tex.c03.IL7R',
                                                              if_else(CD8.Tex.harmony$label == 'Tex.c4','Tex.c04.CREM','aa'))))

# signature score ---------------------------------------------------------

stem.feature = c('CCR7','TCF7','SELL','LEF1','CCR5')
resident.feature = c('NR4A1','NR4A3','CD69','CXCR6','ITGAE')
cyto.feature = c('IFNG','GNLY','GZMB','GZMK','GZMH','GZMA','NKG7','FGFBP2')
exhaust.feature =c('TOX2','SOX4','TIGIT','PDCD1','CTLA4','HAVCR2','LAG3','CXCL13')
costi.feature = c('ICOS','TNFSF14','TNFRSF25','TNFRSF9','CD28','TNFSF4')

sig.feature = c(stem.feature, resident.feature, cyto.feature, exhaust.feature, costi.feature)
sig.list=list(stem = stem.feature,resident = resident.feature,cyto = cyto.feature,exhaust = exhaust.feature,costi = costi.feature)

DoHeatmap(CD8.Tex.harmony,features = sig.feature)


count = CD8.Tex.harmony@assays$RNA@counts
library(GSVA)
gsva_es <- gsva(as.matrix(count), min.sz = 3,sig.list, method=c("ssgsea"), kcdf=c("Poisson")) 
signature_exp<-as.matrix(gsva_es)
saveRDS(signature_exp,file = 'manual.ssgsea.rds')

CD8.Tex.harmony@assays$score$score = signature_exp

saveRDS(CD8.Tex.harmony,file = 'CD8.Tex.harmony.rds')

