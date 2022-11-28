# CD8 paper

library(Seurat)
library(ggplot2)
library(paletteer) 

# figure1 -----------------------------------------------------------------

all = readRDS('all.T.filter.mt.geneblacklist.hsp.minC3F100.harmony.rds')

paletteer_d("khroma::vikO")

