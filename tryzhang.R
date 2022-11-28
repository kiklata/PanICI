# try zhang 2022 nature liver classification
library(Seurat)
library(dplyr)
library(harmony)

ALL.info.kegg <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/obj/ALL.info.kegg.remove.blood.rds")

ALL.info.kegg = ALL.info.kegg %>% NormalizeData() %>% ScaleData(features = rownames(ALL.info.kegg)) %>% 
  FindVariableFeatures() %>% RunPCA() 

ALL.info.kegg$sample.ID = paste0(ALL.info.kegg$sample.timepoint,ALL.info.kegg$sample.ID)

ALL.info.kegg =  RunHarmony(ALL.info.kegg,group.by.vars = c('sample.ID','Study'),assay.use='RNA', plot_convergence = F,theta = c(2.5,1.5), 
                            kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

ALL.info.kegg = ALL.info.kegg %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.2,2,0.2))

ALL.info.kegg = ALL.info.kegg %>% RunTSNE(reduction = "harmony",dims = 1:30) %>% RunUMAP(reduction = "harmony",dims = 1:30)

DimPlot(all.harmony,reduction = 'tsne',group.by = 'RNA_snn_res.1')+NoLegend()
table(all.harmony$RNA_snn_res.1)

# remove small cluster
all.harmony = subset(all.harmony,RNA_snn_res.0.4 %in% c('0','1','2','3'))
Idents(all.harmony) = paste0('C',all.harmony$RNA_snn_res.0.4)
DimPlot(all.harmony,reduction = 'umap')

# findmarker for larger metastate
allmarker = FindAllMarkers(all.harmony)


# combine celltype and metastate 
cluster.n = names(table(all.harmony$scibet.celltype))

CD8.n = cluster.n

cd8.n.n = CD8.n[1]
cd8.m.n = CD8.n[c(2:4,17)]
cd8.em.n = CD8.n[c(5,6)]
cd8.emra.n = CD8.n[7]
cd8.k.n = CD8.n[c(8,9)]
cd8.rm.n = CD8.n[10]
cd8.ex.n = CD8.n[11:14]
cd8.isg.n = CD8.n[15]
cd8.mait.n = CD8.n[16]

all.harmony$scibet.celltype.major = 
          if_else(all.harmony$scibet.celltype %in% cd8.n.n ,'Tn',
          if_else(all.harmony$scibet.celltype %in% cd8.m.n ,'Tm',
          if_else(all.harmony$scibet.celltype %in% cd8.em.n,'Tem',
          if_else(all.harmony$scibet.celltype %in% cd8.emra.n,'Temra',
          if_else(all.harmony$scibet.celltype %in% cd8.k.n,'Tnk',
          if_else(all.harmony$scibet.celltype %in% cd8.rm.n,'Trm',
          if_else(all.harmony$scibet.celltype %in% cd8.ex.n,'Tex',
          if_else(all.harmony$scibet.celltype %in% cd8.isg.n,'Tisg',
          if_else(all.harmony$scibet.celltype %in% cd8.mait.n,'Tmait','na')))))))))
DimPlot(all.harmony,reduction = 'umap')|DimPlot(all.harmony,reduction = 'umap',group.by = 'scibet.celltype.major')

all.harmony$Celltype_Metastate = paste0(all.harmony$RNA_snn_res.1,all.harmony$scibet.celltype.major)

DimPlot(all.harmony,reduction = 'umap',group.by = 'Celltype_Metastate')+NoLegend()

all.harmony$RNA_snn_res.22 = all.harmony$Celltype_Metastate
plot.metastate(all.harmony,res = 22)

gmts = cogena::gmt2list('~/PanCancerICI/CD8TMetabolism/Data/KEGG_metabolism_nc.gmt')
res.list = list()
for (i in 1:85) {
  res = checkGeneSymbols(gmts[[i]])
  if(length(names(table(res$Approved)))>1){
    res.list[[names(gmts)[i]]] = res
  }
}

source("~/step1/rungsea.R", echo=TRUE)
C0 = run.gsea(allmarker,clust = 'C0',
              gmt = '~/PanCancerICI/CD8TMetabolism/Data/KEGG_metabolism_nc.gmt',
              method = 'fgsea')
C1 = run.gsea(allmarker,clust = 'C1',
              gmt = '~/PanCancerICI/CD8TMetabolism/Data/KEGG_metabolism_nc.gmt',
              method = 'fgsea')
C2 = run.gsea(allmarker,clust = 'C2',
              gmt = '~/PanCancerICI/CD8TMetabolism/Data/KEGG_metabolism_nc.gmt',
              method = 'fgsea')
C3 = run.gsea(allmarker,clust = 'C3',
              gmt = '~/PanCancerICI/CD8TMetabolism/Data/KEGG_metabolism_nc.gmt',
              method = 'fgsea')
before = subset(all.harmony,sample.timepoint == 'Before')
plot.metastate(before,res = 22)

meta = before@meta.data

#meta$sample.ID = paste0(meta$sample.timepoint,meta$sample.ID)

ptexpan = as.data.frame(table(meta$sample.ID,meta$scibet.celltype))
ptexpan = tidyr::spread(ptexpan,key = 'Var2',value = 'Freq')
rownames(ptexpan) = ptexpan$Var1
ptexpan$Var1 = NULL
ptexpan$total = apply(ptexpan,1,sum)

for (h in 1:(ncol(ptexpan)-1)) {
  ptexpan[,h] = ptexpan[,h]/ptexpan$total
  
}

left.list = c()

for (i in 1:(ncol(ptexpan)-1)) {
  
  freqs = table(ptexpan[,i] == 0)['FALSE']
  left.list[i] = freqs/nrow(ptexpan)
}

library(ConsensusClusterPlus)

results <- ConsensusClusterPlus(as.matrix(ptexpan[,-ncol(ptexpan)]), maxK = 6,
                                reps = 50, pItem = 0.8,
                                pFeature = 0.8,  
                                clusterAlg = "hc", 
                                distance = "pearson",
                                #title = title,
                                plot = "png")  



#cor.mat = cor(as.matrix(ptexpan[,-ncol(ptexpan)]))

#ComplexHeatmap::pheatmap(cor.mat,cellwidth = 4,cellheight = 4,fontsize = 4)
