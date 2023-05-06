library(Seurat) 
library(patchwork) 


CD8.Tex.harmony <- readRDS("~/CD8.Tex.harmony.rds")
tumor = subset(CD8.Tex.harmony,sample.Tissue == 'Tumor')
Idents(tumor) = tumor$treatment.efficacy
res.markers <- FindAllMarkers(tumor, logfc.threshold = 0) 
Idents(tumor) = tumor$sample.timepoint
time.markers <- FindAllMarkers(tumor, logfc.threshold = 0) 
save(res.markers,time.markers,texmarker.rdata)

res = subset(CD8.Tex.harmony,treatment.efficacy == 'R')
Idents(res) = res$sample.timepoint

object.markers <- FindMarkers(res, ident.1 = 'After',ident.2 = 'Before', logfc.threshold = 0) 
object.markers$names <- rownames(object.markers) 
sig_dge.all <- subset(object.markers, p_val_adj<0.05&abs(avg_log2FC)>0.15) 
#所有差异基因 
#View(sig_dge.all) 


Tex.c1 = subset(res,label == 'Tex.c1')
Tex.c2 = subset(res,label == 'Tex.c2')
Tex.c3 = subset(res,label == 'Tex.c3')

Tex.c1.markers <- FindMarkers(Tex.c1, ident.1 = 'After',ident.2 = 'Before', logfc.threshold = 0) 
Tex.c2.markers <- FindMarkers(Tex.c2, ident.1 = 'After',ident.2 = 'Before', logfc.threshold = 0) 
Tex.c3.markers <- FindMarkers(Tex.c3, ident.1 = 'After',ident.2 = 'Before', logfc.threshold = 0) 


library(dplyr) 
library(ggplot2) 
library(ggrepel) 
library(showtext)
font_add(family = "arial", regular = "software/fonts/arial.ttf")
showtext_auto()
fontuse = 'arial'

regroup = function(obj){
  obj$names <- rownames(obj) 
  
  obj <- obj %>% mutate(Difference = pct.1 - pct.2) 
for (i in 1:nrow(obj)){ 
  if (obj$avg_log2FC[i] >= 0.5 )
    { obj$group[i]='up' } 
  else if(obj$avg_log2FC[i] <= -0.5)
  { obj$group[i]='down' } else { obj$group[i]='no' } } 
  return(obj)
}

plot.point = function(marker){
ggplot(marker, aes(x=100*Difference, y=avg_log2FC)) + 
  geom_point(size=1,aes(color=group)) + 
  scale_color_manual(values=c('blue','grey','red'))+ 
  geom_text_repel(data=subset(marker, group !='no'), aes(label=names), segment.size = 0.25, size=2.5)+ 
  geom_vline(xintercept = 0.0,linetype=2)+ geom_hline(yintercept = 0,linetype=2)+ 
  xlab('Percentage Difference of Cells')+ylab('Log2-fold Change')+
  theme_classic()+NoLegend()+
  theme(text = element_text(family = fontuse))
}

Tex.c1.markers = regroup(Tex.c1.markers)
p1 = plot.point(Tex.c1.markers)

Tex.c2.markers = regroup(Tex.c2.markers)
p2 = plot.point(Tex.c2.markers)

Tex.c3.markers = regroup(Tex.c3.markers)
p3 = plot.point(Tex.c3.markers)

object.markers = regroup(object.markers)
p4 = plot.point(object.markers)

p4|p1|p2|p3

ggsave("TopMarkerVol2.pdf", height=8, width=8)