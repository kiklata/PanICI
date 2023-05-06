library(Seurat)
library(harmony)


selected.all.T <- readRDS("~/PaperCD8/data/selected.all.T.rds")
selected.all.T.meta <- readRDS("~/PaperCD8/data/selected.all.T.meta.rds")
selected.all.T@meta.data = selected.all.T.meta
kegg.gene = readRDS('~/PanCancerICI/CD8TMetabolism/Data/metabolic.gene.name.rds')

CD8 = subset(selected.all.T,manual.celltype.minor == 'CD8.Tex')
CD8.kegg = CD8[kegg.gene,]

CD8.kegg = CD8.kegg %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() 
CD8.kegg =  RunHarmony(CD8.kegg,group.by.vars = 'sample.ID',assay.use='RNA', plot_convergence = F,theta = 4, 
                  kmeans_init_nstart=20, kmeans_init_iter_max=5000) 
CD8.kegg = RunUMAP(CD8.kegg,reduction = 'harmony',dims = 1:30,reduction.name = 'umap_harmony')

CD8.kegg = FindNeighbors(CD8.kegg, reduction = "harmony", dims = 1:30)

CD8.kegg = FindClusters(CD8.kegg,resolution = 0.1)

saveRDS(CD8.kegg,file = 'CD8.Tex.kegg.harmony.rds')