library(Seurat)
library(dplyr)
library(harmony)

library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 2400 * 1024^2)
options(future.rng.onMisuse="ignore")


mycol = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
          "#CCEBC5", "#FFED6F", "#66C2A5", "#FC8D62","#8DA0CB", "#E78AC3")

selected.all.T <- readRDS("~/PaperCD8/data/selected.all.T.rds")
selected.all.T.meta <- readRDS("~/PaperCD8/data/selected.all.T.meta.rds")
selected.all.T@meta.data = selected.all.T.meta

s.gene = cc.genes$s.genes
s.gene[15] = 'CENPU'
g2m.gene = cc.genes$g2m.genes
g2m.gene[16] = 'PIMREG'
g2m.gene[30] = 'JPT1'

selected.all.T <- CellCycleScoring(selected.all.T, s.features = s.gene, g2m.features = g2m.gene) 

CD8 = subset(selected.all.T,manual.celltype.major != 'CD4')

CD8 = CD8 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData(vars.to.regress = c('percent.mt','S.Score','G2M.Score')) %>% RunPCA() 

CD8 =  RunHarmony(CD8,group.by.vars = 'sample.ID',assay.use='RNA', plot_convergence = F,theta = 4, 
                  kmeans_init_nstart=20, kmeans_init_iter_max=5000) 

CD8 = RunUMAP(CD8,reduction = 'harmony',dims = 1:20,reduction.name = 'umap_harmony')

CD8 = FindNeighbors(CD8, reduction = "harmony", dims = 1:30)

CD8 = FindClusters(CD8,method = "igraph",algorithm = 4,resolution = seq(0.2,1,0.1))
saveRDS(CD8,file = 'reCD8.rds')

p1 = DimPlot(reCD8,group.by = 'scibet.celltype.minor')
p2 = DimPlot(reCD8,group.by = 'manual.celltype.minor')

p3 = p1|p2

p.list = list()
for (i in seq(0.2,1,0.1)) {
  res = paste0('RNA_snn_res.',i)
  p.list[[res]] = DimPlot(reCD8,group.by = res)
}

ggsave('oldlabel.png',p3,width = 10,height = 4)

for (i in 1:length(p.list)) {
  ggsave(paste0(i,'.png'),p.list[[i]],width = 6,height = 4)
}

# selected res = 1
reCD8$RNA_snn_res.1 = if_else(reCD8$RNA_snn_res.1 == '21','13',as.character(reCD8$RNA_snn_res.1))#regroup cluster21 = 2 into 13 
reCD8$RNA_snn_res.1 = if_else(reCD8$RNA_snn_res.1 == '17','8',as.character(reCD8$RNA_snn_res.1))#regroup cluster17 = 2509 into 8 
reCD8$RNA_snn_res.1 = if_else(reCD8$RNA_snn_res.1 == '18','4',as.character(reCD8$RNA_snn_res.1))#regroup cluster18 = 2145 into 4 
reCD8$RNA_snn_res.1 = if_else(reCD8$RNA_snn_res.1 == '19','2',as.character(reCD8$RNA_snn_res.1))#regroup cluster19 = 1883 into 2 

FeaturePlot(reCD8,features = 'G2M.Score')
reCD8 = subset(reCD8,RNA_snn_res.1 != '13') # remove S1,G2M cycling cell cluster13

Idents(reCD8) = reCD8$RNA_snn_res.1

saveRDS(reCD8,file = 'reCD8.finished.rds')

p = DimPlot(reCD8,label = T)+NoLegend()
ggsave('res.1.CD8.filter.png',p,width = 6,height = 4) 

res.1.marker = FindAllMarkers(reCD8)
res.1.marker = split(res.1.marker,res.1.marker$cluster)

saveRDS(res.1.marker,file = 'res.1.marker.rds')

reCD8 <- readRDS("~/reCD8.finished.rds")

p1 = DimPlot(reCD8,group.by = 'scibet.celltype.minor')
p2 = DimPlot(reCD8,group.by = 'manual.celltype.minor')
p3 = p1|p2
ggsave('oldlabel.filter.png',p3,width = 10,height = 4)

p4 = FeaturePlot(reCD8,features = c('CCR7'), cols = c('grey90','red3'),raster = F)+NoAxes()+NoLegend()
p5 = FeaturePlot(reCD8,features = c('IL7R'), cols = c('grey90','red3'),raster = F)+NoAxes()+NoLegend()
p6 = FeaturePlot(reCD8,features = c('ZNF683'), cols = c('grey90','red3'),raster = F)+NoAxes()+NoLegend()
p7 = FeaturePlot(reCD8,features = c('PDCD1'), cols = c('grey90','red3'),raster = F)+NoAxes()+NoLegend()
p8 = FeaturePlot(reCD8,features = c('CTLA4'), cols = c('grey90','red3'),raster = F)+NoAxes()+NoLegend()
p9 = FeaturePlot(reCD8,features = c('ISG15'), cols = c('grey90','red3'),raster = F)+NoAxes()+NoLegend()
p10 = FeaturePlot(reCD8,features = c('SLC4A10'), cols = c('grey90','red3'),raster = F)+NoAxes()+NoLegend()
p11 = FeaturePlot(reCD8,features = c('GZMB'), cols = c('grey90','red3'),raster = F)+NoAxes()+NoLegend()
p12 = FeaturePlot(reCD8,features = c('GZMK'), cols = c('grey90','red3'),raster = F)+NoAxes()+NoLegend()
p13 = (p4|p5|p6)/(p7|p8|p9)/(p10|p11|p12)

ggsave('feature.filter.pdf',p13,width = 12,height = 10,device=cairo_pdf)

p14 = FeaturePlot(reCD8,features = c('TRDV2'), cols = c('grey90','red3'),raster = F)+NoAxes()+NoLegend()
p15 = FeaturePlot(reCD8,features = c('TRGV9'), cols = c('grey90','red3'),raster = F)+NoAxes()+NoLegend()
ggsave('gdT.png',p14|p15,width = 12,height = 6)

# annotate

# c1 Tem.c02.GZMK
# c2 Tex.c08.CXCL13+GZMK-
# c3 Tem.c03.GZMK
# c4 Tem.c04.GZMH
# c5 Tm.c05.NR4A1
# c6 Tn.c01.CCR7
# c7 Trm.c07.XCL1.ZNF683
# c8 Tem.c06.PDCD1
# c9 T.c12.NK-like
# c10 T.c09.IFN
# c11 MAIT.c13
# c12 gdT.c14.GNLY
# c14 gdT.c15.TRDV1
# c15 T.c10.STMN1
# c16 gdT.c16.TRDV2
# c20 T.c11.ASPN

reCD8$manual.celltype.minor = if_else(
  reCD8$RNA_snn_res.1 == '1','Tem.c02.GZMK',
  if_else(
    reCD8$RNA_snn_res.1 == '2','Tex.c08.CXCL13',
    if_else(
      reCD8$RNA_snn_res.1 == '3','Tem.c03.HMGB2',
      if_else(
        reCD8$RNA_snn_res.1 == '4','Tem.c04.GZMH',
        if_else(
          reCD8$RNA_snn_res.1 == '5','Tm.c05.NR4A1',
          if_else(
            reCD8$RNA_snn_res.1 == '6','Tn.c01.CCR7',
            if_else(
              reCD8$RNA_snn_res.1 == '7','Trm.c07.ZNF683',
              if_else(
                reCD8$RNA_snn_res.1 == '8','T.c06.MHCII',
                if_else(
                  reCD8$RNA_snn_res.1 == '9','T.c12.NK-like',
                  if_else(
                    reCD8$RNA_snn_res.1 == '10','T.c09.IFN',
                    if_else(
                      reCD8$RNA_snn_res.1 == '11','MAIT.c13',
                      if_else(
                        reCD8$RNA_snn_res.1 == '12','γδT.c14.GNLY',
                        if_else(
                          reCD8$RNA_snn_res.1 == '14','γδT.c15.TRDV1',
                          if_else(
                            reCD8$RNA_snn_res.1 == '15','T.c10.STMN1',
                            if_else(
                              reCD8$RNA_snn_res.1 == '16','γδT.c16.TRDV2','T.c11.ASPN')
                          ))))))))))))))
reCD8$manual.celltype.minor = factor(reCD8$manual.celltype.minor,levels = c('Tn.c01.CCR7','Tem.c02.GZMK','Tem.c03.HMGB2','Tem.c04.GZMH',
                                                                            'Tm.c05.NR4A1','T.c06.MHCII','Trm.c07.ZNF683','Tex.c08.CXCL13',
                                                                            'T.c09.IFN','T.c10.STMN1','T.c11.ASPN','T.c12.NK-like',
                                                                            'MAIT.c13','γδT.c14.GNLY','γδT.c15.TRDV1','γδT.c16.TRDV2'))

reCD8$manual.celltype.major = if_else(
  reCD8$RNA_snn_res.1 %in% c('12','14','16'),'γδT',
  if_else(
    reCD8$RNA_snn_res.1 %in% c('1','3','4'),'Effector memory',
    if_else(
      reCD8$RNA_snn_res.1 %in% c('8'),'MHC II',
    if_else(
      reCD8$RNA_snn_res.1 %in% c('2'),'Exhausted',
      if_else(
        reCD8$RNA_snn_res.1 %in% c('5'),'Memory',
          if_else(
            reCD8$RNA_snn_res.1 %in% c('6'),'Naive',
            if_else(
              reCD8$RNA_snn_res.1 %in% c('15','20'),'Cycling',
              if_else(
                reCD8$RNA_snn_res.1 %in% c('7'),'Resident memory',
                if_else(
                  reCD8$RNA_snn_res.1 %in% c('9'),'NK-like',
                  if_else(
                    reCD8$RNA_snn_res.1 %in% c('10'),'Interferon','MAIT'))))))))))

reCD8$manual.celltype.major = factor(reCD8$manual.celltype.major,levels = c('Naive','Effector memory','Memory','MHC II','Resident memory','Exhausted',
                                                                            'Interferon','Cycling','NK-like','MAIT','γδT'))

p1 = DimPlot(reCD8,group.by = 'manual.celltype.major',label = F,cols = mycol[c(1,2,5,6,7,8,9,10,12,13,14)],raster = F)+NoAxes()
p2 = DimPlot(reCD8,group.by = 'manual.celltype.minor',label = F,cols = mycol,raster = F)+NoAxes()
ggsave('dimplot_manual.pdf', p1|p2,width = 12,height = 4, dpi = 300,device=cairo_pdf)

#p4 = DimPlot(reCD8,group.by = 'manual.celltype.major',label = F,cols = colorqual,split.by = 'sample.timepoint')
#p5 = DimPlot(reCD8,group.by = 'manual.celltype.major',label = F,cols = colorqual,split.by = 'treatment.efficacy')

library(ggplot2)

tumor = subset(reCD8.finished,sample.Tissue == 'Tumor')
blood = subset(reCD8.finished,sample.Tissue == 'Blood')

before = subset(tumor,sample.timepoint == 'Before')
after = subset(tumor,sample.timepoint == 'After')

# celltype.minor by sample before 
source("~/PaperCD8/code/plotprop.R")
p1 = plot.prop(before,cluster = 'manual.celltype.minor',timepoint = 'Pre-Treatment',bar.col = 'turquoise4')
p2 = plot.prop(after,cluster = 'manual.celltype.minor',timepoint = 'Post-Treatment',bar.col = 'turquoise4')
ggsave('prepost_prop.pdf',p1|p2,width = 10,height = 4,device=cairo_pdf)

before.r = subset(before,treatment.efficacy == 'R')
before.nr = subset(before,treatment.efficacy == 'NR')

p3 = plot.prop(before.r,cluster = 'manual.celltype.minor',timepoint = 'Response',bar.col = 'turquoise4')
p4 = plot.prop(before.nr,cluster = 'manual.celltype.minor',timepoint = 'Non Response',bar.col = 'turquoise4')
ggsave('preRNR_prop.pdf',p3|p4,width = 10,height = 4,device=cairo_pdf)

after.r = subset(after,treatment.efficacy == 'R')
after.nr = subset(after,treatment.efficacy == 'NR')

p5 = plot.prop(after.r,cluster = 'manual.celltype.minor',timepoint = 'Response',bar.col = 'turquoise4')
p6 = plot.prop(after.nr,cluster = 'manual.celltype.minor',timepoint = 'Non Response',bar.col = 'turquoise4')
ggsave('postRNR_prop.pdf',p5|p6,width = 10,height = 4,device=cairo_pdf)

# before by study
study.all = names(table(before$Study))[-5]

for (i in 1:length(study.all)) {
  seu.r = subset(before.r,Study == study.all[i])
  seu.nr = subset(before.nr,Study == study.all[i])
  p1 = plot.prop(seu.r,cluster = 'manual.celltype.minor',timepoint = 'Response',bar.col = 'turquoise4')
  p2 = plot.prop(seu.nr,cluster = 'manual.celltype.minor',timepoint = 'Non Response',bar.col = 'turquoise4')
  ggsave(paste0(study.all[i],'.pre.RnR.png'),p1|p2,width = 10,height = 4)
}

# after by study
study.all = names(table(after$Study))

for (i in 1:length(study.all)) {
  seu.r = subset(after.r,Study == study.all[i])
  seu.nr = subset(after.nr,Study == study.all[i])
  p1 = plot.prop(seu.r,cluster = 'manual.celltype.minor',timepoint = 'Response',bar.col = 'turquoise4')
  p2 = plot.prop(seu.nr,cluster = 'manual.celltype.minor',timepoint = 'Non Response',bar.col = 'turquoise4')
  ggsave(paste0(study.all[i],'.post.RnR.png'),p1|p2,width = 10,height = 4)
}

#plot.fc
source("~/PaperCD8/code/plotfc.R")
reCD8 <- readRDS("~/PaperCD8/data/reCD8.finished.rds")

tumor = subset(reCD8,sample.Tissue == 'Tumor')
before = subset(tumor,sample.timepoint == 'Before')

p = plot.fc(before,cluster = 'manual.celltype.minor',norm.col = 'steelblue',sign.col = 'firebrick4',mytitle = 'Pre-Treatment')
ggsave('pretilFC.pdf',p,width = 6,height = 4,device=cairo_pdf)

# bc_bassez----------------
bassez = subset(reCD8,Study == 'BC_Bassez')
bassez$treatment.efficacy = ifelse(bassez$treatment.efficacy=='R','E','NE')
before.ba = subset(bassez,sample.timepoint == 'Before')
after.ba = subset(bassez,sample.timepoint == 'After')
before.ba.e = subset(before.ba,treatment.efficacy == 'E')
before.ba.ne = subset(before.ba,treatment.efficacy == 'NE')
after.ba.e = subset(after.ba,treatment.efficacy == 'E')
after.ba.ne = subset(after.ba,treatment.efficacy == 'NE')

p1 = DimPlot(before.ba.e,group.by = 'manual.celltype.major',
        cols = mycol[c(1,2,5,6,7,8,9,10,12,13,14)],raster = F)+
  NoAxes()+labs(title = '',tag = 'E',x = 'UMAP1',y = 'UMAP2')+
  theme(plot.tag = element_text(face = 'plain',size = 10),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        #axis.title.x = element_text(colour = 'black'),
        axis.title.y = element_text(colour = 'black',angle = 90,size = 10))
p2 = DimPlot(before.ba.ne,group.by = 'manual.celltype.major',
        cols = mycol[c(1,2,5,6,7,8,9,10,12,13,14)],raster = F)+
  NoAxes()+labs(title = '',tag = 'NE')+
  theme(plot.tag = element_text(face = 'plain',size = 10),
        axis.line.x = element_line(size = 0.5),
        #axis.line.y = element_line(size = 0.5),
        #axis.title.x = element_text(colour = 'black')
        )

p3 = p1+p2+plot_layout(guides = 'collect')+
  plot_annotation(title = 'Pre-Treatment',theme = theme(plot.title = element_text(hjust = 0.4)))
p3

p4 = DimPlot(after.ba.e,group.by = 'manual.celltype.major',
             cols = mycol[c(1,2,5,6,7,8,9,10,12,13,14)],raster = F)+
  NoAxes()+labs(title = '',tag = 'E',x = 'UMAP1',y = 'UMAP2')+
  theme(plot.tag = element_text(face = 'plain',size = 10),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        axis.title.x = element_text(colour = 'black',size = 10),
        axis.title.y = element_text(colour = 'black',size = 10,angle = 90))
p5 = DimPlot(after.ba.ne,group.by = 'manual.celltype.major',
             cols = mycol[c(1,2,5,6,7,8,9,10,12,13,14)],raster = F)+
  NoAxes()+labs(title = '',tag = 'NE')+
  theme(plot.tag = element_text(face = 'plain',size = 10),
        axis.line.x = element_line(size = 0.5),
        #axis.line.y = element_line(size = 0.5),
        axis.title.x = element_text(colour = 'black',size = 10))

p6 = p4+p5+plot_layout(guides = 'collect')+
  plot_annotation(title = 'Post-Treatment',theme = theme(plot.title = element_text(hjust = 0.4)))

p6
ggsave('bassez.pre.pdf',p3,width = 10,height = 4,device=cairo_pdf)
ggsave('bassez.post.pdf',p6,width = 10,height = 4,device=cairo_pdf)
