library(Seurat)

setwd('~/PaperCD8/data/Tex')

selected.all.T <- readRDS("~/PaperCD8/data/selected.all.T.rds")

selected.all.T.meta <- readRDS("~/PaperCD8/data/selected.all.T.meta.rds")
selected.all.T@meta.data = selected.all.T.meta

CD8 = subset(selected.all.T,manual.celltype.minor == 'CD8.Tex')

before = subset(CD8,sample.timepoint == 'Before')
after = subset(CD8,sample.timepoint == 'After')

sample = names(table(before$sample.ID))

for(i in 1:length(sample)){
  
  seu = subset(before,sample.ID == sample[i])
  seu <- NormalizeData(seu, normalization.method = "RC",scale.factor = 1e6)
  
  data = as.data.frame(seu@assays$RNA@data)
  write.table(data,file = paste0('compass/before/',sample[i],'.tsv'),
              row.names = T,col.names = T,sep = '\t')
}


sample = names(table(after$sample.ID))

for(i in 1:length(sample)){
  
  seu = subset(after,sample.ID == sample[i])
  seu <- NormalizeData(seu, normalization.method = "RC",scale.factor = 1e6)
  
  data = as.data.frame(seu@assays$RNA@data)
  write.table(data,file = paste0('compass/after/',sample[i],'.tsv'),
              row.names = T,col.names = T,sep = '\t')
}

