# tsne

library(Seurat)
library(dplyr)

ALL.info.kegg <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/ALL.info.kegg.rds")

ALL.info.kegg = FindNeighbors(ALL.info.kegg)
ALL.info.kegg = FindClusters(ALL.info.kegg,resolution = 0.2)


table(ALL.info.kegg$sample.ID)


# pre-----------------------


pre.info.kegg = pre.info.kegg %>% ScaleData(features = rownames(pre.info.kegg)) %>% 
  FindVariableFeatures() %>% RunPCA() 

pre.info.kegg = RunTSNE(pre.info.kegg,reduction.name = 'tsne_pca')

pre.info.kegg = FindNeighbors(pre.info.kegg)
pre.info.kegg = FindClusters(pre.info.kegg,resolution = 0.2)


ptexpan = as.data.frame(table(ALL.info.kegg$sample.ID,ALL.info.kegg$seurat_clusters))
ptexpan = tidyr::spread(ptexpan,key = 'Var2',value = 'Freq')
rownames(ptexpan) = ptexpan$Var1
ptexpan$Var1 = NULL
ptexpan$total = apply(ptexpan,1,sum)

for (i in 1:(ncol(ptexpan)-1)) {
  ptexpan[,i] = ptexpan[,i]/ptexpan$total
}

ptefficacy = as.data.frame(table(ALL.info.kegg$sample.ID,ALL.info.kegg$treatment.efficacy))
ptefficacy = tidyr::spread(ptefficacy,key = 'Var2',value = 'Freq')
rownames(ptefficacy) = ptefficacy$Var1
ptefficacy$Var1 = NULL
ptefficacy$type = if_else(ptefficacy$NR==0,'R','NR')

ptexpan$type = ptefficacy$type

long.pt.pre <- tidyr::gather(ptexpan,key=metastate,value = prop,-c('total','type'))

library(ggplot2)
library(ggsci)
library(ggpubr)

ggboxplot(long.pt.pre, x = "metastate", y = "prop",
          fill = "type", palette = "nejm")+
  stat_compare_means(aes(group = type),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))

