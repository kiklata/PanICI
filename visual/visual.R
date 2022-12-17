
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


# 0.2 func load -----------------------------------------------------------

source("~/PaperCD8/code/plotfc.R")
source("~/PaperCD8/code/clusterbox.R")
source("~/PaperCD8/code/texclusterbox.R")
source("~/PaperCD8/code/majorstack.R")
source("~/PaperCD8/code/plotprop.R")


# 0.3 color set -----------------------------------------------------------

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


# 1.1 CD8 dimplot ---------------------------------------------------------


# 1.2 CD8 featureplot -----------------------------------------------------


# 1.3 CD8 prop plot -------------------------------------------------------


# 1.4 CD8 FC plot ---------------------------------------------------------


# 1.5 CD8 Bassez dimplot --------------------------------------------------


# 1.6 CD8 stackplot -------------------------------------------------------


# 1.7 CD8 avgheatmap ------------------------------------------------------


# 2.1 Tex dimplot ---------------------------------------------------------


# 2.2 Tex dotplot ---------------------------------------------------------


# 2.3 Tex subcluster heatmap ----------------------------------------------


# 2.4 Tex subcluster RvsNR boxplot ----------------------------------------


# 2.5 Tex.c02.GZMH RvsNR boxplot ------------------------------------------


# 2.6 Tex RvsNR forestplot ------------------------------------------------


# 2.7 Tex velcano ---------------------------------------------------------


# 2.8 Tex ssgsea score boxplot -------------------------------------------


# 2.9 Tex avgheatmap ------------------------------------------------------


# 2.10 trajectory ----------------------------------------------------------


# 2.11 Tex scenic ---------------------------------------------------------


# 2.12 Tex compass --------------------------------------------------------


