library(Seurat)
library(dplyr)
setwd('~/PanCancerICI/TumorMetabolism/Data/scData')

# AntiPD1_BC------------------
raw_cohort2_tumor <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/raw_cohort2_tumor.rds")
pre = subset(raw_cohort2_tumor, timepoint == 'Pre')
Idents(pre) = pre$expansion
markers = FindMarkers(pre,ident.1 = 'E',ident.2 = 'NE')
markers = filter(markers,p_val_adj<0.05)
AntiPD1_BC_cohort2_pre_marker = markers

on = subset(raw_cohort2_tumor, timepoint == 'On')
Idents(on) = on$expansion
markers = FindMarkers(on,ident.1 = 'E',ident.2 = 'NE')
markers = filter(markers,p_val_adj<0.05)
AntiPD1_BC_cohort2_on_marker = markers

raw_cohort1_tumor <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/raw_cohort1_tumor.rds")
raw_cohort1_tumor = subset(raw_cohort1_tumor,expansion != 'n/a')
pre = subset(raw_cohort1_tumor, timepoint == 'Pre')
Idents(pre) = pre$expansion
markers = FindMarkers(pre,ident.1 = 'E',ident.2 = 'NE')
markers = filter(markers,p_val_adj<0.05)
AntiPD1_BC_cohort1_pre_marker = markers

on = subset(raw_cohort1_tumor, timepoint == 'On')
Idents(on) = on$expansion
markers = FindMarkers(on,ident.1 = 'E',ident.2 = 'NE')
markers = filter(markers,p_val_adj<0.05)
AntiPD1_BC_cohort1_on_marker = markers
print('AntiPD1_BC')
save(AntiPD1_BC_cohort1_pre_marker,AntiPD1_BC_cohort2_pre_marker,
     AntiPD1_BC_cohort1_on_marker,AntiPD1_BC_cohort2_on_marker,
     file = 'AntiPD1_BC/RvsNR_DEG/Anti_PD1_BC.rdata')

# GSE123813_BCC-----------------------------
raw_BCC_tumor <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/GSE123813_BCC/raw_BCC_tumor.rds")

raw_BCC_tumor$response = if_else(raw_BCC_tumor$patient %in% c('su001','su002','su003','su004','su009','su010'),
                                 'Yes','No')
Idents(raw_BCC_tumor) = raw_BCC_tumor$response
pre = subset(raw_BCC_tumor,treatment == 'pre')
markers = FindMarkers(pre,ident.1 = 'Yes',ident.2 = 'No')
markers = filter(markers,p_val_adj<0.05)
GSE123813_BCC_pre_marker = markers

post = subset(raw_BCC_tumor,treatment == 'post')
markers = FindMarkers(post,ident.1 = 'Yes',ident.2 = 'No')
markers = filter(markers,p_val_adj<0.05)
GSE123813_BCC_post_marker = markers
print('GSE123813_BCC')
save(GSE123813_BCC_pre_marker,GSE123813_BCC_post_marker,file = 'GSE123813_BCC/RvsNR_DEG/GSE123813_BCC.rdata')

# SCP1288_ccRCC------------------------------
raw_tumor <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/SCP1288_ccRCC/raw_tumor.rds")
icb = subset(raw_tumor, ICB_Exposed == 'ICB')
icb$response = if_else(icb$ICB_Response %in% c('ICB_PR','ICB_SD'),'Yes',
                       if_else(icb$ICB_Response == 'ICB_PD','No','NE'))
Idents(icb) = icb$response
markers = FindMarkers(icb,ident.1 = 'Yes',ident.2 = 'No')
markers = filter(markers,p_val_adj<0.05)
SCP1288_ccRCC_marker = markers
print('SCP1288_ccRCC')
save(SCP1288_ccRCC_marker,file = 'SCP1288_ccRCC/RvsNR_DEG/SCP1288_ccRCC.rdata')

