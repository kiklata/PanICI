#Figure 1
library(dplyr)
library(Seurat)

#A-----------------
Before.kegg.final <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/before/Before.kegg.final.rds")

Idents(Before.kegg.final) = Before.kegg.final$RNA_snn_res.0.5

# subset predicted cell, 
Before.kegg.plot = 
  subset(Before.kegg.final,RNA_snn_res.0.5 %in% c('0','1','2','3','4','5','6','7'))


Before.kegg.plot$RNA_snn_res.0.5 = factor(Before.kegg.plot$RNA_snn_res.0.5,levels = c('0','1','2','3','4','5',
                                                                      '6','7'))
Before.kegg.plot = RunTSNE(Before.kegg.plot)
# metastate dimplot
DimPlot(Before.kegg.plot,group.by = 'RNA_snn_res.0.5')
Before.kegg.plot$MetabolicState = Before.kegg.plot$RNA_snn_res.0.5
Idents(Before.kegg.plot) = Before.kegg.plot$MetabolicState

saveRDS(Before.kegg.plot,file = "~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/before/plot/before.kegg.plot.rds")

Before.kegg.plot <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/before/plot/before.kegg.plot.rds")

source("~/PanCancerICI/CD8TMetabolism/Analysis/CD8/plotmetastate.R", echo=TRUE)

Before.kegg.plot = RunTSNE(Before.kegg.plot)

library(SCP)
p1a1 = ClassDimPlot(
  srt = Before.kegg.plot, group.by = c("MetabolicState"),
  reduction = "UMAP", theme_use = "theme_blank",show_stat = F)
p1a1
ggsave(filename = 'p1a1.png',p1a1,width = 6,height = 4)

p1a2 = plot.metastate(obj = Before.kegg.plot,res = 0.5)
p1a2
ggsave(filename = 'p1a2.png',p1a2,width = 6,height = 4)



#subset for cancertype
source("~/PanCancerICI/CD8TMetabolism/Analysis/CD8/plotsuremetastate.R", echo=TRUE)

study.all = names(table(Before.kegg.plot$Study))[-4]
p.list = list()
for (i in 1:length(study.all)) {
  seu = subset(Before.kegg.plot,Study == study.all[i])
  p.list[[study.all[i]]] = plot.sure.metastate(seu,cluster = c('4','5'))
}


# findmarker
all.plotmarker = FindAllMarkers(Before.kegg.plot)
#subset for cluster 4,5,8
all.plotmarker = filter(all.plotmarker,cluster %in% c('4','5'))
all.plotmarker = filter(all.plotmarker,p_val_adj<0.05)

cluster4marker = filter(all.plotmarker,cluster == '4')
cluster5marker = filter(all.plotmarker,cluster == '5')

cluster4marker = filter(cluster4marker,avg_log2FC>0)
cluster5marker = filter(cluster5marker,avg_log2FC>0)

# GO

run.go = function(obj){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  genetoanno = obj$gene
  genelist.ENTREZID = bitr(genetoanno,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,2]

  ego = enrichGO(genelist.ENTREZID,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",
               ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",
               #universe = ,qvalueCutoff = 0.2,
               minGSSize = 10,maxGSSize = 500,readable = FALSE,pool = FALSE)
  return(ego)
}

cluster4go = run.go(cluster4marker)
cluster5go = run.go(cluster5marker)

cluster4cell = colnames(subset(Before.kegg.plot,MetabolicState =='4'))
cluster5cell = colnames(subset(Before.kegg.plot,MetabolicState =='5'))


# kegg


run.kegg = function(obj){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  genetoanno = obj$gene
  genelist.ENTREZID = bitr(genetoanno,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,2]
  
  ekegg = enrichKEGG(genelist.ENTREZID,organism = 'hsa',keyType = "kegg",
                 pvalueCutoff = 0.05,pAdjustMethod = "BH")
  return(ekegg)
}

cluster4kegg = run.kegg(cluster4marker)
cluster5kegg = run.kegg(cluster5marker)

barplot(cluster4kegg)
barplot(cluster5kegg)


cluster4cell = colnames(subset(Before.kegg.plot,MetabolicState =='4'))
cluster5cell = colnames(subset(Before.kegg.plot,MetabolicState =='5'))

# no use
# GSVA
library(scMetabolism)
Before.kegg.plot<-sc.metabolism.Seurat(obj = Before.kegg.plot, method = "GSVA", 
                                        imputation = F, ncores = 1, metabolism.type = "KEGG")
Before.kegg.plot.subset.gsva <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/before/plot/Before.kegg.plot.subset.gsva.rds")

meta.score.gsva = data.frame(cbind(apply(Before.kegg.plot.subset.gsva[,cluster4cell],1,mean),
      apply(Before.kegg.plot.subset.gsva[,cluster5cell],1,mean),
      apply(Before.kegg.plot.subset.gsva[,cluster8cell],1,mean)))
colnames(meta.score.gsva) = c('clust4','clust5','clust8')
# complexheatmap
heatmap(as.matrix(meta.score.gsva),scale = NULL,Colv = NA)
ComplexHeatmap::pheatmap(as.matrix(meta.score.gsva),
                         cluster_rows = F,cluster_cols = F,
                         legend = T, color = mycol,
                         show_colnames = F,angle_col = '0',
                         main = 'GSVA',
                         cellwidth = 5, cellheight = 5,fontsize = 5)

# ssGSEA
library(scMetabolism)
Before.kegg.plot<-sc.metabolism.Seurat(obj = Before.kegg.plot, method = "ssGSEA", 
                                       imputation = F, ncores = 1, metabolism.type = "KEGG")
Before.kegg.plot.subset.ssgsea <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/before/plot/Before.kegg.plot.subset.ssgsea.rds")
Before.kegg.final.all.ssgsea <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/before/Before.kegg.final.all.ssgsea.rds")

meta.score.ssgsea = data.frame(cbind(apply(Before.kegg.plot.subset.ssgsea[,cluster4cell],1,mean),
                                   apply(Before.kegg.plot.subset.ssgsea[,cluster5cell],1,mean),
                                   apply(Before.kegg.plot.subset.ssgsea[,cluster8cell],1,mean)))
colnames(meta.score.ssgsea) = c('clust4','clust5','clust8')

cluster0cell = colnames(subset(Before.kegg.plot,MetabolicState =='0'))
all.meta.score.ssgsea = data.frame(apply(Before.kegg.final.all.ssgsea[,cluster0cell],1,mean))

for (i in 2:length(names(table(Before.kegg.plot$MetabolicState)))) {
  
  cell = colnames(subset(Before.kegg.plot,MetabolicState == paste0(i-1)))
  score = data.frame(apply(Before.kegg.final.all.ssgsea[,cell],1,mean))
  all.meta.score.ssgsea = cbind(all.meta.score.ssgsea,score)
  
}
colnames(all.meta.score.ssgsea) = paste0('cluster',0:9)

heatmap(as.matrix(meta.score.ssgsea),scale = NULL,Colv = NA)

ComplexHeatmap::pheatmap(as.matrix(meta.score.ssgsea[c(carbohydrate,energy,lipid,Aminoacid,glycan),]),
                         cluster_rows = F,cluster_cols = F,
                         legend = T, color = mycol,
                         show_colnames = F,angle_col = '0',
                         main = 'ssGSEA',
                         cellwidth = 5, cellheight = 5,fontsize = 5)



# cor with PanCancer TIL label
mycol = rev(corrplot::COL2('RdBu',n = 200))

ComplexHeatmap::pheatmap(scale(table(Before.kegg.plot$scibet.celltype,Before.kegg.plot$MetabolicState)),
                         cluster_rows = F,cluster_cols = F,
                         legend = F, color = mycol,scale = 'none',
                         show_colnames = T,angle_col = '0',
                         main = 'CellType',#border_color = NA,
                         cellwidth = 10, cellheight = 10,fontsize = 8)

