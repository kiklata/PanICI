# Defined variable
# sample.timepoint,treatment.type,treatment.efficacy,sample.ID
# sample.timepoint: before/after/NonTreat
# treatment.type: PD1/PDL1/CTLA4/PD1_CTLA4/NonTreat
# treatment.efficacy: R/NR/NoInfo based on original study
# sample.ID
# CancerType
# Study
# Platform: 10X/Smart
# DateType: Count/TPM

library(Seurat)
library(dplyr)

# BC_PD1---------
cohort2.on.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.on.T.cd8.rds")
cohort2.pre.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.pre.T.cd8.rds")
cohort1.on.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.on.T.cd8.rds")
cohort1.pre.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.pre.T.cd8.rds")

bc_pd1 = merge(cohort1.pre.T.cd8,c(cohort1.on.T.cd8,cohort2.pre.T.cd8,cohort2.on.T.cd8))

bc_pd1$sample.timepoint = if_else(bc_pd1$timepoint == 'Pre','Before','After')
bc_pd1$treatment.type = 'PD1'
bc_pd1$treatment.efficacy = if_else(bc_pd1$expansion =='E','R','NR')
bc_pd1$sample.ID = bc_pd1$patient_id
saveRDS(bc_pd1,file = "~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/BC_PD1.CD8.rds")

# BC_PDL1------------------

pre.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/pre.T.cd8.rds")
post.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/post.T.cd8.rds")

bc_pdl1 = merge(pre.T.cd8,post.T.cd8)

bc_pdl1$sample.timepoint = if_else(bc_pdl1$Group == 'Pre-treatment','Before','After')
bc_pdl1$treatment.type = 'PDL1'
bc_pdl1$treatment.efficacy = if_else(bc_pdl1$Efficacy == 'PR','R','NR')
bc_pdl1$sample.ID = bc_pdl1$Patient
saveRDS(bc_pdl1,file = "~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/BC_PDL1.CD8.rds")

#BCC_PD1--------------
post.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/post.T.cd8.rds")
pre.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/pre.T.cd8.rds")

bcc_pd1 = merge(pre.T.cd8,post.T.cd8)

bcc_pd1$sample.timepoint = if_else(bcc_pd1$treatment == 'pre','Before','After')
bcc_pd1$treatment.type = 'PD1'
bcc_pd1$treatment.efficacy = if_else(bcc_pd1$patient %in% c('su001','su002', 'su003', 'su004','su009','su012'),'R','NR')
bcc_pd1$sample.ID = bcc_pd1$patient
saveRDS(bcc_pd1,file = "~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/BCC_PD1.CD8.rds")

#ccRCC---------------
post.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/post.T.cd8.rds")
ccrcc = post.T.cd8
ccrcc$sample.timepoint = 'After'
ccrcc$treatment.type = if_else(ccrcc$donor_id %in% c('P915'),'PD1_CTLA4','PD1')
ccrcc$treatment.efficacy = if_else(ccrcc$ICB_Response =='ICB_PR','R','NR')
ccrcc$sample.ID = ccrcc$donor_id
saveRDS(ccrcc,file = "~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/ccRCC.CD8.rds")

#Melanoma_Moshe------------------
pre.ctla <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_CTLA4_Moshe/pre.T.cd8.rds")
pre.dual <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/pre.T.cd8.rds")
post.dual <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/post.T.cd8.rds")
pre.pd1 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/pre.T.cd8.rds")
post.pd1 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/post.T.cd8.rds")

melan = merge(pre.ctla,c(pre.dual,post.dual,pre.pd1,post.pd1))

melan$sample.timepoint = if_else(melan$timepoint == 'Pre','Before','After')
melan$treatment.type = if_else(melan$therapy == 'anti-PD1','PD1',
                               if_else(melan$therapy == 'anti-CTLA4','CTLA4','PD1_CTLA4'))
melan$treatment.efficacy = if_else(melan$Response == 'Responder','R','NR')
melan$sample.ID = melan$PtID
saveRDS(melan,file = "~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/Melanoma_Moshe.CD8.rds")

# Melanoma_livnat---------------
post.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/post.T.cd8.rds")
naive.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/naive.T.cd8.rds")

livnat = merge(post.T.cd8,naive.T.cd8)
livnat$sample.timepoint = if_else(livnat$treatment.group == 'post.treatment','After','NonTreat')
livnat$treatment.type = if_else(livnat$samples %in% c('Mel04.3','Mel58','Mel60','Mel88','Mel98'),'CTLA4',
                                if_else(livnat$samples %in% c('Mel121.1','Mel74'),'PD1',
                                        if_else(livnat$samples %in% c('Mel102','Mel106','Mel110','Mel126','Mel194','Mel72','Mel75',
                                                                      'Mel75.1','Mel78','Mel94'),'PD1_CTLA4','NonTreat')))
livnat$sample.ID = livnat$samples
livnat$treatment.efficacy = 'NoInfo'
# do not have response info
saveRDS(livnat,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/Melanoma_Livnat.CD8.rds')

#nscle_liu--------------

post.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/post.T.cd8.rds")
pre.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/pre.T.cd8.rds")

liu = merge(pre.T.cd8,post.T.cd8)
liu$sample.timepoint = if_else(liu$timepoint == 'Pre','Before','After')
liu$treatment.type = if_else(liu$patient %in% c('P36','P37','P38','P1','P10','P13','P19','P29','P30','P33','P35'),'PD1','NonTreat')
liu$sample.ID = liu$sample
liu$treatment.efficacy = if_else(liu$patient %in% c('P36','P37','P38'),'NR',
                                 if_else(liu$patient %in% c('P1','P10','P13','P19','P29','P30','P33','P35'),'R','NonTreat'))
saveRDS(liu,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/NSCLC_Liu.CD8.rds')

#nscle_caushi_til------------------
T.cd8.list <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Caushi/T.cd8.list.rds")

name.all = character()

for(i in 1:46){
  names = T.cd8.list[[i]]$sample[1]
  name.all[i] = names
}

T.cd8 = merge(x = T.cd8.list[[1]],y = T.cd8.list[2:46],add.cell.ids = name.all)

T.cd8$sample.timepoint = 'After'
T.cd8$treatment.type = 'PD1'
T.cd8$sample.ID = T.cd8$sample
# MPR as R, non-MPR as NR
T.cd8$treatment.efficacy = if_else(substring(T.cd8$sample,1,8) %in% 
                                     c( "MD01-005_tumor_3","MD01-005_tumor_4","MD01-005_tumor_5", "MD01-005_tumor_6",    
                                        "MD01-005_tumor_7","MD01-005_tumor_8","MD01-005_tumor_9","MD01-010_tumor_1",
                                        "MD043-003_tumor_1","MD043-003_tumor_2","MD043-003_tumor_3","MD043-003_tumor_4",   
                                        "MD043-003_tumor_5","MD043-008_tumor_1","NY016-022_tumor_1","NY016-022_tumor_2",   
                                        "NY016-022_tumor_3","NY016-022_tumor_4","NY016-025_tumor_1","NY016-025_tumor_2",
                                        "NY016-025_tumor_3","NY016-025_tumor_4"),'R','NR')
saveRDS(T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Caushi/NSCLC.cd8.rds')

# scc---------
post.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/post.T.cd8.rds")
pre.T.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/pre.T.cd8.rds")

scc = merge(pre.T.cd8,post.T.cd8)

scc$sample.timepoint = if_else(scc$treatment == 'pre','Before','After')
scc$treatment.type = 'PD1'
scc$treatment.efficacy = if_else(scc$patient %in% c('su010','su011'),'R','NR')
scc$sample.ID = scc$patient
saveRDS(scc,file = "~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/SCC_PD1.CD8.rds")

# merge all dataset-------------------

BC_PD1.CD8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/BC_PD1.CD8.rds")
BC_PDL1.CD8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/BC_PDL1.CD8.rds")
BCC_PD1.CD8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/BCC_PD1.CD8.rds")
ccRCC.CD8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/ccRCC.CD8.rds")
Melanoma_Moshe.CD8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/Melanoma_Moshe.CD8.rds")
Melanoma_Livnat.CD8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/Melanoma_Livnat.CD8.rds")
NSCLC.cd8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Caushi/NSCLC.cd8.rds")
NSCLC_Liu.CD8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/NSCLC_Liu.CD8.rds")
SCC_PD1.CD8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/SCC_PD1.CD8.rds")

BC_PD1.CD8$CancerType = 'Breast'
BC_PDL1.CD8$CancerType = 'Breast'
BCC_PD1.CD8$CancerType = 'BCC'
ccRCC.CD8$CancerType = 'ccRCC'
Melanoma_Moshe.CD8$CancerType = 'Melanoma'
Melanoma_Livnat.CD8$CancerType = 'Melanoma'
NSCLC.cd8$CancerType = 'NSCLC'
NSCLC_Liu.CD8$CancerType = 'NSCLC'
SCC_PD1.CD8$CancerType = 'SCC'

BC_PD1.CD8$Study = 'BC_Bassez'
BC_PDL1.CD8$Study = 'BC_Zhang'
BCC_PD1.CD8$Study = 'BCC_Yost'
ccRCC.CD8$Study = 'RCC_Bi'
Melanoma_Moshe.CD8$Study = 'Melanoma_Feldman'
Melanoma_Livnat.CD8$Study = 'Melanoma_Arnon'
NSCLC.cd8$Study = 'NSCLC_Caushi'
NSCLC_Liu.CD8$Study = 'NSCLC_Liu'
SCC_PD1.CD8$Study = 'SCC_Yost'

BC_PD1.CD8$Platform = '10X'
BC_PDL1.CD8$Platform = '10X'
BCC_PD1.CD8$Platform = '10X'
ccRCC.CD8$Platform = '10X'
Melanoma_Moshe.CD8$Platform = 'SmartV2'
Melanoma_Livnat.CD8$Platform = 'SmartV2'
NSCLC.cd8$Platform = '10X'
NSCLC_Liu.CD8$Platform = '10X'
SCC_PD1.CD8$Platform = '10X'

BC_PD1.CD8$DateType = 'Count'
BC_PDL1.CD8$DateType = 'Count'
BCC_PD1.CD8$DateType = 'Count'
ccRCC.CD8$DateType = 'Count'
Melanoma_Moshe.CD8$DateType = 'TPM'
Melanoma_Livnat.CD8$DateType = 'Count'
NSCLC.cd8$DateType = 'Count'
NSCLC_Liu.CD8$DateType = 'Count'
SCC_PD1.CD8$DateType = 'Count'

dataset.name = c('BC_PD1','BC_PDL1','BCC_PD1','ccRCC','Melanoma_Moshe','Melanoma_Livnat','NSCLC_Caushi','NSCLC_Liu','SCC_PD1')
ALL = merge(BC_PD1.CD8,y = c(BC_PDL1.CD8,BCC_PD1.CD8,ccRCC.CD8,Melanoma_Moshe.CD8,
                                 Melanoma_Livnat.CD8,NSCLC.cd8,NSCLC_Liu.CD8,SCC_PD1.CD8),add.cell.ids = dataset.name)

ALL$sample.ID = paste0(ALL$Study,ALL$sample.ID)
meta = ALL@meta.data
meta = meta[,c(2,3,11,13:20)]
ALL@meta.data = meta
saveRDS(ALL,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/All.rds')

