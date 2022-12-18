
# 0.1 pkg load ----------------------------------------------------------

library(Seurat)
library(dplyr)

library(ggplot2)
library(ggplotify)
library(ggthemes)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(ggsci)
library(patchwork)
library(showtext)
font_add('arial',regular = '~/software/fonts/arial.ttf')
showtext_auto()

library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 2400 * 1024^2)
options(future.rng.onMisuse="ignore")

# 0.2 color set -----------------------------------------------------------

# Naive effector.memory memory MHCII resi.memory exhausted interferon cycling nk.like mait gdt
CD8.dim.major.col = c('#8dd3c7','#ffffb3','#80b1d3','#8da0cb','#b3de69','#fb8072','#bebada','#d9d9d9','#fc8d62','#ccebc5','#fccde5')
CD8.dim.minor.col = c('#8dd3c7','#ffffb3','#fdb462','#ffed6f','#80b1d3','#8da0cb','#b3de69','#fb8072','#bebada','#d9d9d9','#66c2a5','#fc8d62','#ccebc5','#fccde5','#bc80bd','#e78ac3')

prop.col = '#00868b'

FC.high.col = '#8b1a1a'
FC.low.col = '#4682b4'
FC.line.col = '#bebebe'

CD8.stack.col = CD8.dim.major.col

heatmap.col = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7",
                                 "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)

Tex.dim.col = c('#e41a1c','#377eb8','#4daf4a','#984ea3')

# bassez zhang liu caushi rcc bcc scc luoma
study.col = c('#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6')

Tex.bar.col = c('#4dbbd5','#3c5488','#f39b7f','#00a087')

RvsNR.box.col = c('#bf2c2d','#1b91bb')

# E.pre NE.pre E.post NE.post
Texc02.box.col = c('#de7597','#bf2c2d','#4e8abc','#1b91bb')

# R.post R.pre NR.post NR.pre 
Tex.volcano.point.col = c('#ff3030','#009acd','#cd853f','#2e8b57')
Tex.volcano.margin.col = c('#4d4d4d')

Tex.score.box.col = Tex.dim.col


# 0.3 obj load ------------------------------------------------------------
CD8 <- readRDS("~/PaperCD8/data/reCD8.finished.rds")

Tprop <- read.csv("~/PaperCD8/data/Tprop.csv")

CD8.Tex.harmony <- readRDS("~/PaperCD8/data/Tex/CD8.Tex.harmony.rds")
Tex.scenic.rss <- readRDS("~/PaperCD8/data/Tex/Tex.scenic.rss.rds")
volcanodf <- readRDS("~/PaperCD8/data/Tex/volcanodf.rds")

# 1.1 CD8 dimplot ---------------------------------------------------------


# 1.2 CD8 featureplot -----------------------------------------------------
gene = c('CCR7','IL7R','ZNF683','CXCL13','PDCD1','ISG15','SLC4A10','GZMB','GZMK')

# 1.3 CD8 prop plot -------------------------------------------------------
source("~/PaperCD8/code/plotprop.R")


# 1.4 CD8 FC plot ---------------------------------------------------------
source("~/PaperCD8/code/plotfc.R")


# 1.5 CD8 Bassez dimplot --------------------------------------------------


# 1.6 CD8 stackplot -------------------------------------------------------
source("~/PaperCD8/code/majorstack.R")


# 1.7 CD8 avgheatmap ------------------------------------------------------
gene = c('CCR7','SELL','LEF1','TCF7','IL7R',
         'GZMK','ITM2C','CST7','CD28','KLRG1',
         'NR4A1','EGR1','BAG3','UBE2S',
         'HLA-DRA','HLA-DRB1','HLA-DQA1','HLA-DRB5','HLA-DPA1',
         'XCL1','ZNF683','ANXA1','VIM','LMVA',
         'CXCL13','KRT86','LAYN','GZMB','CTLA4','PDCD1','HAVCR2','TIGIT',
         'IFI6','MX1','ISG15','IFIT3','IFIT1',
         'STMN1','TUBB','LINC00861','S100A4','HMGN2',
         'FGFBP2','KKLRD1','KLRC3','KLRF1','TYROBP',
         'KLRB1','CCL20','CEBPD','CCR6','SLC4A10',
         'TRDC','TRDV2','TRGV9')

# 2.1 Tex dimplot ---------------------------------------------------------


# 2.2 Tex dotplot ---------------------------------------------------------
gene = c('NR4A2','ZNF331','NR4A3','CXCR4','MYADM','CREM', # c04
         'TNFRSF4','KLRB1','SELL','CCR7','LMNA','LTB','IL7R', # c03
         'KLRG1','FCGR3A','CX3CR1','GZMH','NKG7','PLEK','FGFBP2', #c02
         'CXCL13','TNFRSF9','LAG3','IFNG','GZMB','CCL4','CCL4L2','CCL3' #c01 
         )

# 2.3 Tex subcluster heatmap ----------------------------------------------
gene = c('NR4A2','ZNF331','NR4A3','CXCR4','MYADM','CREM', # c04
         'TNFRSF4','KLRB1','SELL','CCR7','LMNA','LTB','IL7R', # c03
         'HAVCR2','CD63','CXCR6','CX3CR1','ZNF683','GNLY','GZMH','CCL5', #c02
         'CXCL13','TNFRSF9','LAG3','IFNG','GZMB','CCL4','CCL4L2','CCL3' #c01
         )


# 2.4 Tex subcluster RvsNR boxplot ----------------------------------------
source("~/PaperCD8/code/clusterbox.R")


# 2.5 Tex.c02.GZMH RvsNR boxplot ------------------------------------------
source("~/PaperCD8/code/texclusterbox.R")


# 2.6 Tex RvsNR forestplot ------------------------------------------------


# 2.7 Tex velcano ---------------------------------------------------------


# 2.8 Tex ssgsea score boxplot -------------------------------------------


# 2.9 Tex avgheatmap ------------------------------------------------------
stem.feature = c('CCR7','TCF7','SELL','LEF1','CCR5')
resident.feature = c('NR4A1','NR4A3','CD69','CXCR6','ITGAE')
cyto.feature = c('IFNG','GNLY','GZMB','GZMK','GZMH','GZMA','NKG7','FGFBP2')
exhaust.feature =c('TOX2','SOX4','TIGIT','PDCD1','CTLA4','HAVCR2','LAG3','CXCL13')
costi.feature = c('ICOS','TNFSF14','TNFRSF25','TNFRSF9','CD28','TNFSF4')

gene = c(stem.feature, resident.feature, cyto.feature, exhaust.feature, costi.feature)

# 2.10 trajectory ----------------------------------------------------------


# 2.11 Tex scenic ---------------------------------------------------------


# 2.12 Tex compass --------------------------------------------------------


