library(Seurat)
library(dplyr)

score.hngc <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/score.hngc.rds")
ALL.info.kegg.remove.blood <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/obj/ALL.info.kegg.remove.blood.rds")
before = subset(ALL.info.kegg.remove.blood,sample.timepoint == 'Before')

before.ssgsea = score.hngc[,colnames(score.hngc) %in% colnames(before)]

cor.mat1 = cor(before.ssgsea,method = 'spearman')
hc1 = hclust(as.dist(1 - cor.mat1), method = 'ward.D2')


saveRDS(hc1,file= 'before.hc.rds')

before.hc <- readRDS("~/before.hc.rds")
k.select =4
tree = as.data.frame(cutree(before.hc,k = k.select))
colnames(tree)[1] = 'k.select'



table(tree$k.select)

#Before.kegg.final <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/before/Before.kegg.final.rds")
Before.kegg.final = AddMetaData(Before.kegg.final,tree)
#Before.kegg.final = RunTSNE(Before.kegg.final)

source("~/PanCancerICI/CD8TMetabolism/Analysis/CD8/plotmetastate.R", echo=TRUE)

Before.kegg.final$RNA_snn_res.2 = Before.kegg.final$k.select
plot.metastate(obj = Before.kegg.final,res = 2)

Before.kegg.final$RNA_snn_res.3 = Before.kegg.final$scibet.celltype
plot.metastate(obj = Before.kegg.final,res = 3)




DimPlot(Before.kegg.final,group.by = 'k.select')
mycol = rev(corrplot::COL2('RdBu',n = 200))

hccelltype = as.data.frame(table(Before.kegg.final$scibet.celltype,Before.kegg.final$k.select))
hccelltype = tidyr::spread(hccelltype,key = 'Var2',value = 'Freq')
rownames(hccelltype) = hccelltype$Var1
hccelltype$Var1 = NULL
hccelltype$total = apply(hccelltype,1,sum)
hccelltype$all = sum(hccelltype$total)

for (h in 1:(ncol(hccelltype)-2)) {
  hccelltype[,h] = (hccelltype[,h]/hccelltype$total)/(hccelltype$total/hccelltype$all)
}

ComplexHeatmap::pheatmap(hccelltype[,c(-5,-6)],
                         cluster_rows = F,cluster_cols = F,
                         legend = F, color = mycol,scale = 'none',
                         show_colnames = T,angle_col = '0',
                         main = 'CellType',#border_color = NA,
                         cellwidth = 10, cellheight = 10,fontsize = 8)




source("~/PanCancerICI/CD8TMetabolism/Analysis/CD8/plotsuremetastate.R", echo=TRUE)

study.all = names(table(Before.kegg.final$Study))
p.list = list()
for (i in 1:length(study.all)) {
  seu = subset(Before.kegg.plot,Study == study.all[i])
  p.list[[study.all[i]]] = plot.sure.metastate(seu,cluster = c('4','5'))
}


