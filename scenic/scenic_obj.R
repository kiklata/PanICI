library(dplyr)
library(Seurat)
library(SeuratDisk)
library(SCENIC)

setwd('~/PaperCD8/data/Tex/scenic')

seu = readRDS('~/CD8.Tex.harmony.rds')

count = seu@assays$RNA@counts
nCountsPerGene = apply(count,MARGIN = 1,sum)
count[count>0] = 1
nCellsPerGene = apply(count,MARGIN = 1,sum)
  
nCells= ncol(seu)
minCountsPerGene = round(.05*nCells,2) # 0.05 counts per cell
minSamples=round(.05*nCells,0) # 5% of cells
  
GeneSelectCounts = nCountsPerGene[nCountsPerGene>minCountsPerGene]
GeneSelectSamples = nCellsPerGene[nCellsPerGene>minSamples]
GeneSelect = intersect(names(GeneSelectCounts),names(GeneSelectSamples))
  
exp = seu@assays$RNA@counts
exp = exp[rownames(exp) %in% GeneSelect,]
seufilter = CreateSeuratObject(exp)
seufilter = AddMetaData(seufilter,metadata = seu@meta.data[,c(-1,-2,-3)])
SaveH5Seurat(seufilter, filename = "filter.h5Seurat",assay = 'RNA')
Convert("filter.h5Seurat", dest = "h5ad")
  

