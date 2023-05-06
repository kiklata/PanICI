library(Seurat)
library(dplyr)
source("~/PanCancerICI/CD8TMetabolism/Analysis/CD8/plotmetastate.R", echo=TRUE)

ALL.info.kegg <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/ALL.info.kegg.remove.blood.rds")

all.time = names(table(ALL.info.kegg$sample.timepoint))

all.time.plot = list()
for (i in 1:length(all.time)) {
  
  time.info.kegg = subset(ALL.info.kegg,sample.timepoint == all.time[i])
  
  time.info.kegg = time.info.kegg %>% ScaleData(features = rownames(time.info.kegg)) %>% 
    FindVariableFeatures() %>% RunPCA() 
  
  time.info.kegg = RunTSNE(time.info.kegg)
  
  time.info.kegg = FindNeighbors(time.info.kegg)
  time.info.kegg = FindClusters(time.info.kegg,resolution = seq(0.1,1,0.1))
  
  plot.list = list()
  
  for (k in 1:10) {
    res = as.character(k/10)
    p1 = plot.metastate(obj = time.info.kegg,res = res)
    plot.list[[res]] = p1
  }
  
  all.time.plot[[all.time[i]]] = plot.list
  saveRDS(time.info.kegg,file = paste0('~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/',all.time[i],'.tsne.kegg.rds'))
}

saveRDS(all.time.plot,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/all.time.plot.rds')
