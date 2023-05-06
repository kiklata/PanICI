library(Seurat)
ALL.info.hgnc <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/obj/ALL.info.hgnc.rds")
source("~/step1/geneblacklist.R", echo=TRUE)

rpgene = grep('^RP[LS]',rownames(ALL.info.hgnc),value = T)
final.gene =setdiff(rownames(ALL.info.hgnc),y = c(immgene,tcrgene,rpgene))
ALL.info.hgnc = ALL.info.hgnc[final.gene,]

library(SeuratDisk)
SaveH5Seurat(ALL.info.hgnc, filename = "ALL.info.hgnc.h5Seurat",assay = 'RNA')
Convert("ALL.info.hgnc.h5Seurat", dest = "h5ad")
