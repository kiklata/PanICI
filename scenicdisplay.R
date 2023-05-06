library(SCENIC)


CD8.Tex.harmony <- readRDS("~/CD8.Tex.harmony.rds")

meta = CD8.Tex.harmony@meta.data

aucres = as.data.frame(read_csv('~/PaperCD8/data/Tex/scenic/res/Tex.aucell.csv'))
rownames(aucres) = aucres$Cell
aucres$Cell = NULL
CD8.Tex.harmony = AddMetaData(CD8.Tex.harmony,aucres)

colnames(CD8.Tex.harmony@meta.data)[33:136] = gsub(pattern = "\\...",replacement = "(+)",colnames(CD8.Tex.harmony@meta.data)[33:136])

FeaturePlot(CD8.Tex.harmony,features = 'FOS(+)',cols = c('grey90','red3'))
VlnPlot(CD8.Tex.harmony,features = 'FOS(+)',group.by = 'label',pt.size = 0)

aucres = as.matrix(t(aucres))
rss = as.data.frame(calcRSS(AUC = aucres, cellAnnotation = CD8.Tex.harmony$label))
rss=na.omit(rss) 

fc = apply(rss, 1, sum)-apply(rss, 1, median)  

topgene = names(fc[order(fc,decreasing = T)][1:20])

select.gene = c('AR(+)','CEBPB(+)','ATF3(+)','FOXO1(+)','HES1(+)','HMGA1(+)','RB1(+)','NR2F2(+)',
                'SOX4(+)','SOX9(+)','TCF3(+)','STAT1(+)','BRCA1(+)','E2F1(+)','JUN(+)','NFKB1(+)')

Stand = function(data){
  new.data = (data-min(data))/(max(data)-min(data))
  #new.data = if_else(new.data > 0.5,1,0)
}
plot.rss = rss[select.gene,] %>% Stand() %>% as.matrix()

mycol = rev(corrplot::COL2('RdBu',200))

ComplexHeatmap::pheatmap(plot.rss,
                         cluster_rows = F,cluster_cols = F,
                         legend = T,
                         color = mycol,
                         show_colnames = T,angle_col = '0',
                         main = 'Tumor')
