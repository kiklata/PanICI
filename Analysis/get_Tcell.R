# BC_PD1_cohort1-----------------------
raw_cohort1 <- readRDS("~/PanCancerICI/Data/scDATA/AntiPD1_BC/raw_cohort1.rds")
raw_cohort1 = subset(raw_cohort1,expansion !='n/a')
pre = subset(raw_cohort1,timepoint == 'Pre')
on = subset(raw_cohort1,timepoint =='On')
pre.t = subset(pre,cellType == 'T_cell')
on.t = subset(on,cellType == 'T_cell')
saveRDS(pre.t,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.pre.T.rds')
saveRDS(on.t,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.on.T.rds')
# BC_PD1_cohort2-------------------------
raw_cohort2 <- readRDS("~/PanCancerICI/Data/scDATA/AntiPD1_BC/raw_cohort2.rds")
raw_cohort2 = subset(raw_cohort2,expansion !='n/a')
pre2 = subset(raw_cohort2,timepoint == 'Pre')
on2 = subset(raw_cohort2,timepoint =='On')
pre2.t = subset(pre2,cellType == 'T_cell')
on2.t = subset(on2,cellType == 'T_cell')
saveRDS(pre2.t,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.pre.T.rds')
saveRDS(on2.t,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.on.T.rds')

# BC_PDL1-----------------------------
scBC <- readRDS("~/PanCancerICI/Data/scDATA/AntiPDL1_BC_immune/RNA/scBC.rds")
scBC = subset(scBC,Treatment == 'Anti-PD-L1+Chemo')
scBC = subset(scBC,Major.celltype == 'T cell')
pre = subset(scBC,Group == 'Pre-treatment')
post = subset(scBC,Group == 'Post-treatment')
saveRDS(pre,'~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/pre.T.rds')
saveRDS(post,'~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/post.T.rds')

# NSCLC_Caushi------------------------
# all TIL

# NSCLC_Liu----------------
# all TIL
raw <- readRDS("~/PanCancerICI/Data/scDATA/Liu_NSCLC_TIL/raw.rds")
raw$timepoint = if_else(raw$sample %in% 
                   c("P1.post.1","P1.post.2", "P1.post.3","P10.post.1",  
 "P13.post.1", "P13.post.2","P19.post.1","P29.post.1","P30.post.1" ,
 "P33.post.1","P35.post.1", "P36.post.1" ,"P37.post.1" ,"P38.post.1"),'Post','Pre')

pre = subset(raw,timepoint =='Pre')
post = subset(raw,timepoint =='Post')

saveRDS(pre,'~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/pre.T.rds')
saveRDS(post,'~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/post.T.rds')

# Melanoma_Livnat--------------------------

raw = subset(raw,cell.types %in% c('T.CD4','T.CD8','T.cell'))

naive = subset(raw,treatment.group =='treatment.naive')
post = subset(raw,treatment.group =='post.treatment')

saveRDS(naive,'~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/naive.T.rds')
saveRDS(post,'~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/post.T.rds')

# Melanoma_Moshe-----------------------
# all immune cell
# need annotate
tpm <- readRDS("~/PanCancerICI/Data/scDATA/GSE120575_Melanoma_immune/tpm.rds")

tpm$timepoint = if_else(tpm$PtID %in% 
                          c("Post_P1","Post_P1_2","Post_P11", "Post_P12","Post_P13",
                            "Post_P14","Post_P15","Post_P16", "Post_P17", "Post_P18",
                            "Post_P19","Post_P2","Post_P20","Post_P21","Post_P22",
                            "Post_P23","Post_P23_2","Post_P28","Post_P28_2","Post_P3",
                            "Post_P3_2", "Post_P4", "Post_P6","Post_P8"),'Post','Pre')
pre = subset(tpm,timepoint == 'Pre')
post = subset(tpm,timepoint == 'Post')
pre.pd1 = subset(pre,therapy == 'anti-PD1')
pre.ctla = subset(pre,therapy == 'anti-CTLA4')
pre.dual = subset(pre,therapy == 'anti-CTLA4+PD1')
post.pd1 = subset(post,therapy == 'anti-PD1')
post.dual = subset(post,therapy == 'anti-CTLA4+PD1')

saveRDS(pre.pd1,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/pre.rds')
saveRDS(post.pd1,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/post.rds')
saveRDS(pre.ctla,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_CTLA4_Moshe/pre.rds')
saveRDS(pre.dual,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/pre.rds')
saveRDS(post.dual,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/post.rds')

# BCC---------------
raw_BCC <- readRDS("~/PanCancerICI/Data/scDATA/GSE123813_BCC_SCC/raw_BCC.rds")
raw_BCC = subset(raw_BCC,cluster %in% c('CD8_mem_T_cells','CD4_T_cells','CD8_act_T_cells',
                                        'CD8_ex_T_cells','Tcell_prolif','Tregs'))

pre = subset(raw_BCC, treatment == 'pre')
post = subset(raw_BCC, treatment == 'post')

saveRDS(pre,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/pre.T.rds')
saveRDS(post,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/post.T.rds')

# SCC------------------
raw_SCC <- readRDS("~/PanCancerICI/Data/scDATA/GSE123813_BCC_SCC/raw_SCC.rds")
pre = subset(raw_SCC,treatment == 'pre')
post = subset(raw_SCC,treatment == 'post')
saveRDS(pre,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/pre.T.rds')
saveRDS(post,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/post.T.rds')

# ccRCC_Kevin--------------

raw <- readRDS("~/PanCancerICI/Data/scDATA/SCP1288_ccRCC/raw.rds")

raw = subset(raw,FinalCellType %in% c('41BB-Hi CD8+ T cell','41BB-Lo CD8+ T cell',
                                      'Cycling CD8+ T cell','MitoHigh CD8+ T cell',
                                      'MX1-Hi CD8+ T cell','Effector T-Helper',
                                      'Memory T-Helper','MitoHigh T-Helper','NKT',
                                      'T-Reg'))

raw = subset(raw,ICB_Exposed == 'ICB')
raw = subset(raw,ICB_Response != 'ICB_NE')
saveRDS(raw,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/post.rds')
