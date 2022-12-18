library(Seurat)
library(dplyr)
obj=readRDS("~/PaperCD8/data/reCD8.finished.rds")
CD8.Tex.harmony <- readRDS("~/CD8.Tex.harmony.rds")

CD8.Tex.label = as.data.frame(CD8.Tex.harmony$label)
c1cell = rownames(filter(CD8.Tex.label,`CD8.Tex.harmony$label` == 'Tex.c1'))
c2cell = rownames(filter(CD8.Tex.label,`CD8.Tex.harmony$label` == 'Tex.c2'))
c3cell = rownames(filter(CD8.Tex.label,`CD8.Tex.harmony$label` == 'Tex.c3'))
c4cell = rownames(filter(CD8.Tex.label,`CD8.Tex.harmony$label` == 'Tex.c4'))

CD8 = subset(obj,manual.celltype.minor %in% c('Tn.c01.CCR7','Tem.c02.GZMK','Tm.c05.NR4A1','Tem.c04.GZMH',
                                              'T.c06.MHCII','Trm.c07.ZNF683','Tex.c08.CXCL13'))

CD8$manual.celltype.minor = if_else(colnames(CD8) %in% c1cell,'Tex.c1',
                                    if_else(colnames(CD8) %in% c2cell,'Tex.c2',
                                            if_else(colnames(CD8) %in% c3cell,'Tex.c3',
                                                    if_else(colnames(CD8) %in% c4cell,'Tex.c4',
                                                            as.character(CD8$manual.celltype.minor)))))


# downsample -------------------------------------------------------------------
Idents(CD8) = CD8$manual.celltype.minor
CD8.downsample = subset(CD8,downsample = 1000)
# each celltype minor 1000
#CD8.downsample = CD8
