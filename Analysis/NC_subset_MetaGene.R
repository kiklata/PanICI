library(dplyr)
library(Seurat)

metabolic.gene.name <- readRDS("~/PanCancerICI/TumorMetabolism/Data/metabolic.gene.name.rds")

# BC-------------------
raw_cohort1_tumor <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/raw_cohort1_tumor.rds")
raw_cohort1_tumor = raw_cohort1_tumor[metabolic.gene.name,]

raw_cohort1_tumor = subset(raw_cohort1_tumor,expansion !='n/a')
pre = subset(raw_cohort1_tumor,timepoint == 'Pre')
on = subset(raw_cohort1_tumor,timepoint =='On')

score.cohort1.hc.pre <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/score.cohort1.hc.pre.rds")
score.cohort1.hc.on <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/score.cohort1.hc.on.rds")

load("~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/score.cohort1.rdata")

# PRE
pre = FindVariableFeatures(pre) %>% NormalizeData() %>% ScaleData(features = rownames(pre)) %>% 
  RunPCA() %>% RunTSNE(dims = 1:20) %>% FindNeighbors() %>% FindClusters(resolution = 0.02)

DimPlot(pre,reduction = 'tsne')

k.hclust = 6
tree = as.data.frame(cutree(score.cohort1.hc.pre,k = k.hclust))
colnames(tree) = 'hcres'
table(tree$hcres)

pre = AddMetaData(pre,tree)
pre$seurat_clusters = if_else(pre$seurat_clusters == '0','TSNE_c1',
                              if_else(pre$seurat_clusters == '1','TSNE_c2',
                                      if_else(pre$seurat_clusters == '2','TSNE_c3',
                                              if_else(pre$seurat_clusters == '3','TSNE_c4',
                                                      if_else(pre$seurat_clusters == '4','TSNE_c5','TSNE_c6'
                                                              )))))
pre$hcres = if_else(pre$hcres == '1','HC_c1',
                              if_else(pre$hcres == '2','HC_c2',
                                      if_else(pre$hcres == '3','HC_c3',
                                              if_else(pre$hcres == '4','HC_c4',
                                                      if_else(pre$hcres == '5','HC_c5','HC_c6'
                                                      )))))

heatmap = pheatmap::pheatmap(scale(table(pre$hcres,pre$seurat_clusters)),cluster_rows = F,cluster_cols = F)

pre$seurat_clusters = factor(pre$seurat_clusters,levels = c('TSNE_c1','TSNE_c3','TSNE_c5',
                                                            'TSNE_c2','TSNE_c4','TSNE_c6'))
heatmap = pheatmap::pheatmap(scale(table(pre$hcres,pre$seurat_clusters)),cluster_rows = F,cluster_cols = F)

library(tidyr)
hcpred = data.frame(table(pre$expansion,pre$hcres))
hcpred = spread(hcpred,key = 'Var2',value = 'Freq')
rownames(hcpred) = hcpred$Var1
hcpred$Var1 = NULL
hcpred$sum = apply(hcpred, 1, sum)

for (i in 1:(ncol(hcpred)-1)) {
  hcpred[,i] = hcpred[,i]/hcpred$sum
}

hcpred[3,] = hcpred[1,]-hcpred[2,]

# find hc_c1,hc_c2 as key cell state
DimPlot(pre,cells.highlight = hc_c2_cell)

hc_c1_cell = rownames(filter(tree,hcres == '1'))
nohc1 = rownames(filter(tree,hcres != '1'))
hc_c2_cell = rownames(filter(tree,hcres == '2'))

score.cohort1.hc.pre.hc1 = data.frame(score.cohort1.pre.norm[hc_c1_cell,])
score.cohort1.hc.pre.hc2 = data.frame(score.cohort1.pre.norm[hc_c2_cell,])
score.cohort1.hc.pre.hc1$clust = 'c1'
score.cohort1.hc.pre.hc2$clust = 'c2'

score.cohort1.hc.pre.nohc1 = data.frame(score.cohort1.pre.norm[nohc1,])
score.cohort1.hc.pre.nohc1$clust = 'noc1'

#long.score.pre <- gather(rbind(score.cohort1.hc.pre.hc1,score.cohort1.hc.pre.nohc1)
#                         ,key=ssGSEA,value = Expression,-c('clust'))
#library(ggsci)
#library(ggpubr)

#ggboxplot(long.score.pre, x = "ssGSEA", y = "Expression",
#          fill = "clust", palette = "lancet")+
#  stat_compare_means(aes(group = clust),
#                     method = "wilcox.test",
#                     label = "p.signif",
#                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
#                                      symbols = c("***", "**", "*", "ns")))

res.table = list()
for (i in 1:(ncol(score.cohort1.hc.pre.hc1)-1)) {
  res = t.test(x = score.cohort1.hc.pre.hc1[,i],y = score.cohort1.hc.pre.nohc1[,i])
  res.a = data.frame(p.val = res$p.value,hc1.mean = res[["estimate"]][["mean of x"]],nohc1.mean =res[["estimate"]][["mean of y"]])
  rownames(res.a) = colnames(score.cohort1.hc.pre.hc1)[i]
  res.table[[colnames(score.cohort1.hc.pre.hc1)[i]]] = res.a
  
}

res.hc1 = data.table::rbindlist(res.table)
rownames(res.hc1) = names(res.table)
res.hc1$p.val.adj = p.adjust(res.hc1$p.val,method = 'BH')
res.hc1 = filter(res.hc1,p.val.adj<0.01)
res.hc1$pathway = rownames(res.hc1)
res.hc1$fc = res.hc1$hc1.mean - res.hc1$nohc1.mean



# by pts

long.pt.pre <- tidyr::gather(hcpts[,-7],key=cellstate,value = prop,-c('type'))

library(tidyr)
hcpts = data.frame(table(pre$patient_id,pre$hcres))
hcpts = spread(hcpts,key = 'Var2',value = 'Freq')
rownames(hcpts) = hcpts$Var1
hcpts$Var1 = NULL
hcpts$sum = apply(hcpts, 1, sum)

for (i in 1:(ncol(hcpts)-1)) {
  hcpts[,i] = hcpts[,i]/hcpts$sum
}

ptexpan = as.data.frame(table(pre$patient_id,pre$expansion))
ptexpan = spread(ptexpan,key = 'Var2',value = 'Freq')
ptexpan$type = if_else(ptexpan$E>0,'E','NE')
rownames(ptexpan) = ptexpan$Var1
ptexpan$Var1 = NULL

hcpts$type = ptexpan$type

ggboxplot(hcpts,x = 'type',y = 'HC_c1')

ggboxplot(long.pt.pre, x = "cellstate", y = "prop",
                 fill = "type", palette = "lancet")+
           stat_compare_means(aes(group = type),
                               method = "wilcox.test",
                               label = "p.signif",
                               symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                symbols = c("***", "**", "*", "ns")))

logit.res = list()
for (i in 1:(ncol(hcpts)-2)) {
  res = summary(glm(as.factor(type)~hcpts[,i],family = binomial(link='logit'),data = hcpts))
  res = res[["coefficients"]]
  logit.res[[colnames(hcpts)[i]]] = res
}
