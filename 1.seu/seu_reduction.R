library(Seurat)
library(dplyr)
library(harmony)


all <- readRDS("~/all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony.rds")

all = all %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData(vars.to.regress = 'percent.mt' ) %>% RunPCA() 

all =  RunHarmony(all,group.by.vars = 'sample.ID',assay.use='RNA', plot_convergence = F,theta = 4, 
                            kmeans_init_nstart=20, kmeans_init_iter_max=5000) 

all = RunUMAP(all,reduction = 'harmony',dims = 1:30,reduction.name = 'umap_harmony')
#all = RunUMAP(all,reduction = 'pca',dims = 1:30,reduction.name = 'umap_pca')

all = FindNeighbors(all, reduction = "harmony", dims = 1:30)
#all = FindNeighbors(all, reduction = "pca", dims = 1:30, graph.name = c('RNA_nn_pca','RNA_snn_pca'))

all = FindClusters(all,method = "igraph",algorithm = 4,resolution = 2)
#harmony.res = all$seurat_clusters

#all = FindClusters(all,method = "igraph",algorithm = 4,resolution = 1,graph.name = 'RNA_snn_pca')
#pca.res = all$seurat_clusters

saveRDS(all,file = 'all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony.rds')


# kegg --------------------------------------------------------------------


library(Seurat)
library(dplyr)
library(harmony)


all <- readRDS("PaperCD8/data/all.T.filter.mt.geneblacklist.hsp.minC3F100.kegg.rds")

all = all %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() 
all =  RunHarmony(all,group.by.vars = 'sample.ID',assay.use='RNA', plot_convergence = F,theta = 4, 
                  kmeans_init_nstart=20, kmeans_init_iter_max=5000) 
all = RunUMAP(all,reduction = 'harmony',dims = 1:30,reduction.name = 'umap_harmony')
#all = RunUMAP(all,reduction = 'pca',dims = 1:30,reduction.name = 'umap_pca')

all = FindNeighbors(all, reduction = "harmony", dims = 1:30)
#all = FindNeighbors(all, reduction = "pca", dims = 1:30, graph.name = c('RNA_nn_pca','RNA_snn_pca'))

all = FindClusters(all,method = "igraph",algorithm = 4,resolution = seq(0.2,1,0.2))
#harmony.res = all$seurat_clusters

saveRDS(all,file = 'all.T.filter.mt.geneblacklist.hsp.minC3F100.kegg.harmony.rds')
