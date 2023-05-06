library(Seurat)
library(dplyr)
library(harmony)
source("~/PanCancerICI/CD8TMetabolism/Analysis/CD8/plotmetastate.R", echo=TRUE)

# maybe two round use metacell first
# raw.batch ---------------------------------------------------------------

Before.tsne.kegg <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/obj/Before.tsne.kegg.rds")
Before.tsne.kegg = subset(Before.tsne.kegg,Study != 'BC_Bassez')

Before.tsne.kegg = Before.tsne.kegg %>% ScaleData(features = rownames(Before.tsne.kegg)) %>% 
  FindVariableFeatures() %>% RunPCA() 

Before.tsne.kegg =  RunHarmony(Before.tsne.kegg,group.by.vars = c('sample.ID','Study'),assay.use='RNA', plot_convergence = F,theta = c(2.5,1.5), 
                        kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

Before.tsne.kegg = Before.tsne.kegg %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.1,1,0.1)) %>% RunTSNE(reduction = "harmony")

Before.tsne.kegg = RunUMAP(Before.tsne.kegg,features = rownames(Before.tsne.kegg))
Before.tsne.kegg = FindClusters(Before.tsne.kegg)

plot.list = list()

for (k in 1:10) {
  res = as.character(k/10)
  p1 = plot.metastate(obj = Before.tsne.kegg,res = res)
  plot.list[[res]] = p1
}

study.all = names(table(Before.tsne.kegg$Study))[-3]
Before.tsne.kegg$MetabolicState = Before.tsne.kegg$RNA_snn_res.0.5
p.list = list()
for (i in 1:length(study.all)) {
  seu = subset(Before.tsne.kegg,Study == study.all[i])
  p.list[[study.all[i]]] = plot.sure.metastate(seu,cluster = c('0'))
}
ComplexHeatmap::pheatmap(scale(table(Before.tsne.kegg$scibet.celltype,Before.tsne.kegg$MetabolicState)),
                         cluster_rows = F,cluster_cols = F,
                         legend = F, color = mycol,scale = 'none',
                         show_colnames = T,angle_col = '0',
                         main = 'CellType',#border_color = NA,
                         cellwidth = 10, cellheight = 10,fontsize = 8)


After.tsne.kegg <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/After.tsne.kegg.rds")

#Before.tsne.kegg = Before.tsne.kegg %>% ScaleData(features = rownames(Before.tsne.kegg)) %>% 
#  FindVariableFeatures() %>% RunPCA() 

After.tsne.kegg =  RunHarmony(After.tsne.kegg,group.by.vars = c('sample.ID','Study'),assay.use='RNA', plot_convergence = F,theta = c(2.5,1.5), 
                               kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

After.tsne.kegg = After.tsne.kegg %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.1,1,0.1)) %>% RunTSNE(reduction = "harmony")
plot.list = list()

for (k in 1:10) {
  res = as.character(k/10)
  p1 = plot.metastate(obj = After.tsne.kegg,res = res)
  plot.list[[res]] = p1
}
