library(Seurat)
library(dplyr)

raw_cohort1 <- readRDS("~/PanCancerICI/Data/scDATA/AntiPD1_BC/raw_cohort1.rds")
load("~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/AntiPD1_BC_pre_NC_ssGSEA_calc.rdata")

raw_cohort1 = subset(raw_cohort1,expansion !='n/a')
all.pre = raw_cohort1[,colnames(pre)]
all.pre = AddMetaData(all.pre,pre@meta.data)

#all.pre$ifhc1 = if_else(all.pre$hcres == 'HC_c1','HC_c1','NO')
#Idents(all.pre) = all.pre$ifhc1

#meta.marker = FindAllMarkers(all.pre)

#saveRDS(meta.marker,file = '~/PanCancerICI/Data/scDATA/AntiPD1_BC/hc_c1.marker.rds')

#print('hc1')

all.pre$ifhc2 = if_else(all.pre$hcres == 'HC_c2','HC_c2','NO')
Idents(all.pre) = all.pre$ifhc2
meta.marker = FindAllMarkers(all.pre)
saveRDS(meta.marker,file = '~/PanCancerICI/Data/scDATA/AntiPD1_BC/hc_c2.marker.rds')
