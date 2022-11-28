library(Seurat)
library(dplyr)
library(SCP)
library(ggplot2)
library(RColorBrewer)

# color setting -----------------------------------------------------------


colors <- subset(brewer.pal.info,category=='qual')
colorqual <- c()
for(i in nrow(colors):1) {
  colorqual <- c(colorqual, brewer.pal(colors[i,'maxcolors'], name = rownames(colors)[i]) )}

# All T cell --------------------------------------------------------------


all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony <- readRDS("~/all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony.rds")

p1 = DimPlot(all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony,label = T,
             cols = colorqual)
p1
ggsave('dimplot_seuratclus.pdf', p1,width = 8,height = 4)

p2 = DimPlot(all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony,group.by = 'scibet.celltype',label = F,cols = colorqual)
p2
ggsave('dimplot_scibet.pdf', p2,width = 12,height = 4)

t.feature = c('CD8A','NKG7','CD4','TRDC','CCR7','CXCL13')
p3 = FeaturePlot(all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony,features = t.feature)
p3
ggsave('feature.pdf',p3,width = 6,height = 8)

p4 = DotPlot(all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony,features = t.feature)
ggsave('dotplot.pdf',p4,width = 4,height = 8)

AverageExpression(all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony,
                  features = t.feature,slot = 'data')

meta = all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony@meta.data
meta$celltype.major = if_else(meta$seurat_clusters %in% 
                                c(1,3,6,7,8,10,12,15,17,20,22,24,26,27,29,31,35,39,40,42,44,46),'CD8',
                              if_else(meta$seurat_clusters %in% 
                                        c(2,4,5,9,11,13,14,16,18,19,25,32,33,36,37,38,41,43,45,47),'CD4',
                                      if_else(meta$seurat_clusters %in% 
                                                c(21,23),'γδ',
                                              if_else(meta$seurat_clusters %in% c(28,30,34),'Mix','NotDefined'))))

# c 28,30,34 might be mix CD4/CD8 NME1 cell proliferated

meta$cluster = if_else(meta$seurat_clusters %in% c('39','40'),'1',
                       if_else(meta$seurat_clusters %in% c('41','43'),'2',
                               if_else(meta$seurat_clusters %in% c('42'),'22',
                                       if_else(meta$seurat_clusters %in% c('43','44','45'),'4',
                                               if_else(meta$seurat_clusters %in% c('46','47'),'10',
                                                       as.character(meta$seurat_clusters))))))
saveRDS(meta,file = 'celltype.major.rds')
all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony@meta.data = meta
Idents(all) = meta$cluster

allmarker = FindAllMarkers(all)
marker.list = split(allmarker, allmarker$cluster)

saveRDS(marker.list,file = 'all.marker.list.rds')

readRDS('allmarker.rds')
all <- readRDS("~/PaperCD8/data/all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony.rds")
celltype.major <- readRDS("~/PaperCD8/data/celltype.major.rds")
all = AddMetaData(all,celltype.major)

p5 = DimPlot(all,group.by = 'celltype.major',label = T,
             cols = colorqual)
ggsave('dimplot_celltypemajor.png', p5,width = 6,height = 4)

# scibet predict --------------------------------------------------------------

library(scibet)
SciBet.ref <- readRDS("~/software/scibet/count/SciBet.ref.rds")

sample.num = names(table(all$sample.ID))

prd.list = list()
for (i in 1:length(sample.num)) {
  
  seu = subset(all,sample.ID == sample.num[i])
  
  expr = seu[['RNA']]@counts
  expr = scuttle::calculateTPM(expr)
  expr = t(as.matrix(expr))
  
  prd = SciBet(SciBet.ref, expr)
  prd = as.data.frame(prd)
  prd$cellid = rownames(expr)
  
  colnames(prd) = c('scibet.celltype','cell.id')
  rownames(prd) = prd$cell.id
  prd.list[[sample.num[i]]] = prd
  print(formattable::percent(i/length(sample.num)))
  
}
saveRDS(prd.list,file = 'prd.list.rds')

prd.table = data.table::rbindlist(prd.list)
celltype.major = celltype.major[,c(1:11,13,14,15,16,17,18)]
rownames(prd.table) = prd.table$cell.id
celltype.major$cell.id = rownames(celltype.major)
celltype.major = left_join(celltype.major,prd.table,by = 'cell.id')
rownames(celltype.major) = celltype.major$cell.id

celltype.major$scibet.celltype.major = substring(celltype.major$scibet.celltype,1,3)

scibet_group <- read_excel("scibet.group.xlsx")
all.minor = names(table(scibet_group$minor))

for (i in 1:length(all.minor)) {
  celltype.major$scibet.celltype.minor = if_else(celltype.major$scibet.celltype %in% 
                                                   c(filter(scibet_group,minor == all.minor[i])[,1])[[1]],all.minor[i],celltype.major$scibet.celltype.minor)
}



# mix ---------------------------------------------------------------------
mix = subset(all,celltype.major == 'Mix')

library(scibet)
SciBet.ref <- readRDS("~/software/scibet/count/SciBet.ref.rds")

expr = mix[['RNA']]@counts
expr = scuttle::calculateTPM(expr)
expr = t(as.matrix(expr))

prd <- SciBet(SciBet.ref, expr)
prd = as.data.frame(prd)
prd$cellid = rownames(expr)

colnames(prd) = c('scibet.celltype.subset','cell.id')
rownames(prd) = prd$cell.id
prd$cell.id = NULL

mix = AddMetaData(mix,prd)

mix$celltype.major = if_else(substring(mix$scibet.celltype.subset,1,3) == 'CD4','CD4','CD8')

DimPlot(mix,group.by = 'celltype.major')|FeaturePlot(mix,features = 'CD8A')

# regroup MIX to CD4/CD8
meta= celltype.major
mix.meta = mix@meta.data
meta = filter(meta,celltype.major !='Mix')
mix.meta = mix.meta[,c(1:18)]
meta.all = rbind(meta,mix.meta)

all = AddMetaData(all,meta.all)
saveRDS(all@meta.data,file = 'celltype.major.rds')

all.scibet.minor = names(table(all$scibet.celltype.minor))
p.list = list()
for (i in 1:length(all.scibet.minor)) {
  p.list[[all.scibet.minor[i]]] = DimPlot(all,cells.highlight = WhichCells(all,idents = all.scibet.minor[i]))
}

for (i in 1:length(p.list)) {
  ggsave(paste0(names(p.list)[i],'.pdf'),p.list[[i]],width = 6,height = 4)
}


# manula_cluster ----------------------------------------------------------
celltype.major$cluster = as.character(celltype.major$cluster)

celltype.major$manual.celltype.minor = 
  if_else(celltype.major$cluster %in% c('28','30','34'),'Cycling',
  if_else(celltype.major$cluster %in% c('27','35','20'),'Deleted',
  if_else(celltype.major$cluster %in% c('21','23','17'),'γδ',
  if_else(celltype.major$cluster %in% c('24'),'CD8.ISG',
  if_else(celltype.major$cluster %in% c('25'),'CD4.ISG',
  if_else(celltype.major$cluster %in% c('26'),'CD8.MAIT',
  if_else(celltype.major$cluster %in% c('12','29','6'),'CD8.Tm',
  if_else(celltype.major$cluster %in% c('1','8','22'),'CD8.Tem',
  if_else(celltype.major$cluster %in% c('7','31'),'CD8.Tex',
  if_else(celltype.major$cluster %in% c('19','16','18'),'CD4.Tfh',
  if_else(celltype.major$cluster %in% c('11','36'),'CD4.Tm',
  if_else(celltype.major$cluster %in% c('9','4'),'CD4.Tem',
  if_else(celltype.major$cluster %in% c('2','10','38'),'CD4.Tn',
  if_else(celltype.major$cluster %in% c('15'),'CD8.Tn',
  if_else(celltype.major$cluster %in% c('5','13','37','14','32','33'),'CD4.Treg',
  if_else(celltype.major$cluster %in% c('3'),'CD8.Trm','NOT'))))))))))))))))

celltype.major$manual.celltype.major = substring(celltype.major$manual.celltype.minor,1,3)
saveRDS(celltype.major,file = "~/PaperCD8/data/celltype.major.rds")
all@meta.data = celltype.major

celltype.major$manual.celltype.major = if_else(celltype.major$manual.celltype.major == 'Cyc',celltype.major$celltype.major,celltype.major$manual.celltype.major)
celltype.major$manual.celltype.minor = if_else(celltype.major$manual.celltype.minor == 'Cycling',paste0(celltype.major$manual.celltype.major,'.Proliferating'),celltype.major$manual.celltype.minor)

p2 = DimPlot(all,group.by = 'manual.celltype.major',label = F,cols = colorqual)

# selected_obj final
all = subset(all,manual.celltype.major != 'Del')

all$manual.celltype.major = if_else(all$manual.celltype.major == 'Cyc',all$celltype.major,all$manual.celltype.major)
all$manual.celltype.minor = if_else(all$manual.celltype.minor == 'Cycling',paste0(all$manual.celltype.major,'.Proliferating'),all$manual.celltype.minor)

saveRDS(all,file = '~/PaperCD8/data/selected.all.T.rds')

p1 = DimPlot(all,group.by = 'manual.celltype.major',label = F,cols = colorqual)
p2 = DimPlot(all,group.by = 'manual.celltype.minor',label = F,cols = colorqual)
p3 = DimPlot(all,group.by = 'cluster',label = T,cols = colorqual)+NoLegend()
ggsave('dimplot_manual.png', p3/(p1|p2),width = 12,height = 10)

p4 = DimPlot(all,group.by = 'manual.celltype.minor',label = F,cols = colorqual,split.by = 'treatment.efficacy')
ggsave('dimplot_RvsNR.png', p4, width = 10,height = 4)


# annotate celltype finished


# annotate metastate ------------------------------------------------------










source("~/step1/rungsea.R", echo=TRUE)

C1 = run.gsea(allmarker,clust = '10',
              gmt = '~/PanCancerICI/CD8TMetabolism/Data/KEGG_metabolism_nc.gmt',
              method = 'fgsea')

marker.list = split(allmarker, allmarker$cluster)