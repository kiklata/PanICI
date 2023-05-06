# find markers

library(Seurat)
library(dplyr)

ALL.info.kegg <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/ALL.info.kegg.remove.blood.rds")

Idents(ALL.info.kegg) = ALL.info.kegg$treatment.efficacy

before = subset(ALL.info.kegg,sample.timepoint == 'Before')
after = subset(ALL.info.kegg,sample.timepoint == 'After')

marker.pre = FindMarkers(before,ident.1 = 'R','NR')
marker.post = FindMarkers(after,ident.1 = 'R','NR')

save(marker.pre,marker.post,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/kegg.markers.rdata')
