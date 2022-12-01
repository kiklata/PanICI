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

CD8.Tex.harmony = FindNeighbors(CD8.Tex.harmony, reduction = "harmony", dims = 1:6)

CD8.Tex.harmony = FindClusters(CD8.Tex.harmony,resolution = seq(0.1,1,0.1))

Idents(CD8.Tex.harmony) = CD8.Tex.harmony$RNA_snn_res.0.1

CD8.Tex.harmony$RNA_snn_res.2 = NULL
CD8.Tex.harmony@meta.data[,c(25:33)] = NULL
CD8.Tex.harmony$seurat_clusters = CD8.Tex.harmony$RNA_snn_res.0.1

saveRDS(CD8.Tex.harmony,file = 'CD8.Tex.harmony.rds')

CD8.Tex.harmony$label = ifelse(CD8.Tex.harmony$RNA_snn_res.0.1 == '0','Tex.c1',
                               if_else(CD8.Tex.harmony$RNA_snn_res.0.1 == '1','Tex.c2','Tex.c3'))
p = DimPlot(CD8.Tex.harmony,raster = F,label = F,group.by = 'label',
        cols = c( "#E41A1C", "#377EB8", "#4DAF4A"))+
  NoAxes()+labs(title = 'Tex',x = 'UMAP1',y = 'UMAP2')+
  theme(plot.title = element_text(face = 'plain',size = 15,hjust = 0.5),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        axis.title.x = element_text(colour = 'black',size = 10),
        axis.title.y = element_text(colour = 'black',angle = 90,size = 10))
ggsave('tex.pdf',p,width = 6,height = 4,device=cairo_pdf)


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

# gene avg.mat
library(pheatmap)
avg.mat = as.matrix(AverageExpression(CD8.Tex.harmony,features = sig.feature)[[1]])
colnames(avg.mat) = c('Tex.c1','Tex.c2','Tex.c3')
pdf('avg.heatmap.pdf',width = 3,height = 8,bg = 'white')
pheatmap(avg.mat,
        cluster_rows = F,cluster_cols = F,
        cellwidth = 15,cellheight = 15,fontsize = 10,
        scale = 'row',border = F,legend_breaks = c(-1,0,1),angle_col = 45,legend = T,
        #annotation_names_row = F,annotation_row = myannorow,
        #gaps_row = c(6,10,18,24),annotation_legend = F,annotation_colors = list(type = c('Stem' = 'grey40',
        #                                                                        'Resident' = 'grey40',
        #                                                                        'Cyto' = 'grey40',
        #                                                                        'Exhausted' = 'grey40',
        #                                                                        'Co-Sti' = 'grey40')),
        )
dev.off()

# ssgsea score
score = CD8.Tex.harmony@assays$score$score
score = as.data.frame(t(score))
CD8.Tex.harmony = AddMetaData(CD8.Tex.harmony,score)

plot.score = CD8.Tex.harmony@meta.data[,c("stem","resident","cyto","exhaust","costi",'label')]
library(ggpubr)
plot.score$label = factor(plot.score$label,levels = c('Tex.c1','Tex.c2','Tex.c3'))
my_comparisons <- list( c("Tex.c2", "Tex.c1"), c("Tex.c3", "Tex.c1"))

mytheme = theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5),
        axis.ticks.y = element_blank())
p.stem = ggboxplot(data = plot.score,x = 'label',y = 'stem',color = 'label',width = 0.5,
                   palette = c("#E41A1C", "#377EB8", "#4DAF4A"),outlier.shape = NA,add = 'none')+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+NoLegend()+
  labs(x = '',y ='Stemness')+mytheme

p.resident = ggboxplot(data = plot.score,x = 'label',y = 'resident',color = 'label',width = 0.5,
                   palette = c("#E41A1C", "#377EB8", "#4DAF4A"),outlier.shape = NA,add = 'none')+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+NoLegend()+
  labs(x = '',y ='Resident')+mytheme
p.cyto = ggboxplot(data = plot.score,x = 'label',y = 'cyto',color = 'label',width = 0.5,
                   palette = c("#E41A1C", "#377EB8", "#4DAF4A"),outlier.shape = NA,add = 'none')+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+NoLegend()+
  labs(x = '',y ='Cytotoxicity')+mytheme
p.exhaust = ggboxplot(data = plot.score,x = 'label',y = 'exhaust',color = 'label',width = 0.5,
                   palette = c("#E41A1C", "#377EB8", "#4DAF4A"),outlier.shape = NA,add = 'none')+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+NoLegend()+
  labs(x = '',y ='Exhaustion')+mytheme
p.costi = ggboxplot(data = plot.score,x = 'label',y = 'costi',color = 'label',width = 0.5,
                      palette = c("#E41A1C", "#377EB8", "#4DAF4A"),outlier.shape = NA,add = 'none')+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+NoLegend()+
  labs(x = '',y ='Co-Stimulatory')+mytheme
p2 = p.stem|p.resident|p.cyto|p.exhaust|p.costi
ggsave('score.box.pdf',p2,width = 10,height = 4,device=cairo_pdf)

# forest plot
