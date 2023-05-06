# Supplementary Table 1

ALL.info.kegg <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/ALL.info.kegg.remove.blood.rds")

meta.all = ALL.info.kegg@meta.data

rownames(meta.all) = NULL

meta.sample = meta.all[,4:11][!duplicated(meta.all$sample.ID),]
write.csv(meta.sample,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/meta.sample.csv',row.names = F)

pts.timepoint = tidyr::spread(as.data.frame(table(ALL.info.kegg$sample.ID,ALL.info.kegg$sample.timepoint)),key = 'Var2',value = 'Freq')

write.csv(pts.timepoint,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/tmp_res/pt.timepoint.csv',row.names = F)
