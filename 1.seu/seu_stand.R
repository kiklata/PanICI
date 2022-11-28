# use all T cell try zhang

# bc_bassez ---------------------------------------------------------------

cohort1.on.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.on.T.rds")
cohort1.pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.pre.T.rds")
cohort2.on.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.on.T.rds")
cohort2.pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.pre.T.rds")

bc_bass = merge(cohort1.pre.T,c(cohort1.on.T,cohort2.pre.T,cohort2.on.T))

bc_bass$sample.timepoint = if_else(bc_bass$timepoint == 'Pre','Before','After')
bc_bass$treatment.type = 'PD1'
bc_bass$treatment.efficacy = if_else(bc_bass$expansion =='E','R','NR')
bc_bass$CancerType = 'BC'
bc_bass$Study = 'BC_Bassez'
bc_bass$Platform = '10X'
bc_bass$DateType = 'Count'
bc_bass$sample.Tissue = 'Tumor'
bc_bass$sample.ID = paste0(bc_bass$Study,bc_bass$sample.timepoint,bc_bass$patient_id)

bc_bass@meta.data[,c(1,4:10)] = NULL

cohort1.on.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.on.T.scibet.rds")
cohort1.pre.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.pre.T.scibet.rds")
cohort2.on.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.on.T.scibet.rds")
cohort2.pre.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.pre.T.scibet.rds")

scibet = rbind(cohort1.on.T.scibet,cohort1.pre.T.scibet,cohort2.on.T.scibet,cohort2.pre.T.scibet)
scibet$cell.id = NULL

bc_bass = AddMetaData(bc_bass,metadata = scibet)

source("~/genecorrect.R", echo=TRUE)
bc_ba = gene.correct(bc_bass)

saveRDS(bc_ba,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/BC_Bassez.rds')

#bc_zhang--------

pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/pre.T.rds")
post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/post.T.rds")


bc_zhang = merge(pre.T,post.T)

bc_zhang$sample.timepoint = if_else(bc_zhang$Group == 'Pre-treatment','Before','After')
bc_zhang$treatment.type = 'PDL1'
bc_zhang$treatment.efficacy = if_else(bc_zhang$Efficacy == 'PR','R','NR')
bc_zhang$CancerType = 'BC'
bc_zhang$Study = 'BC_Zhang'
bc_zhang$Platform = '10X'
bc_zhang$sample.Tissue = if_else(bc_zhang$Origin == 'b','Blood','Tumor')
bc_zhang$DateType = 'Count'
bc_zhang$sample.ID = paste0(bc_zhang$Study,bc_zhang$sample.timepoint,bc_zhang$Sample)

bc_zhang@meta.data[,c(1,4:14)] = NULL

pre.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/pre.T.scibet.rds")
post.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/post.T.scibet.rds")
scibet = rbind(pre.T.scibet,post.T.scibet)

scibet$cell.id = NULL

bc_zhang = AddMetaData(bc_zhang,metadata = scibet)
source("~/genecorrect.R", echo=TRUE)
bc_zh = gene.correct(bc_zhang)

saveRDS(bc_zh,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/BC_Zhang.rds')

# BCC_pd1--------------
post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/post.T.rds")
pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/pre.T.rds")

bcc_pd1 = merge(pre.T,post.T)

bcc_pd1$sample.timepoint = if_else(bcc_pd1$treatment == 'pre','Before','After')
bcc_pd1$treatment.type = 'PD1'
bcc_pd1$treatment.efficacy = if_else(bcc_pd1$patient %in% c('su001','su002', 'su003', 'su004','su009','su012'),'R','NR')
bcc_pd1$CancerType = 'BCC'
bcc_pd1$Study = 'BCC_Yost'
bcc_pd1$Platform = '10X'
bcc_pd1$sample.Tissue = 'Tumor'
bcc_pd1$DateType = 'Count'
bcc_pd1$sample.ID = paste0(bcc_pd1$Study,bcc_pd1$sample.timepoint,bcc_pd1$patient)

bcc_pd1@meta.data[,c(1,4:10)] = NULL

post.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/post.T.scibet.rds")
pre.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/pre.T.scibet.rds")
scibet = rbind(pre.T.scibet,post.T.scibet)
scibet$cell.id = NULL

bcc_pd1 = AddMetaData(bcc_pd1,metadata = scibet)

source("~/genecorrect.R", echo=TRUE)
bcc = gene.correct(bcc_pd1)

saveRDS(bcc,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/BCC.rds')


#SCC_PD1-----------


post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/post.T.rds")
pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/pre.T.rds")

scc = merge(pre.T,post.T)

scc$sample.timepoint = if_else(scc$treatment == 'pre','Before','After')
scc$treatment.type = 'PD1'
scc$treatment.efficacy = if_else(scc$patient %in% c('su010','su011'),'R','NR')
scc$CancerType = 'SCC'
scc$Study = 'SCC_Yost'
scc$Platform = '10X'
scc$sample.Tissue = 'Tumor'
scc$DateType = 'Count'
scc$sample.ID = paste0(scc$Study,scc$sample.timepoint,scc$patient)

scc@meta.data[,c(1,4:10)] = NULL

post.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/post.T.scibet.rds")
pre.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/pre.T.scibet.rds")
scibet = rbind(pre.T.scibet,post.T.scibet)
scibet$cell.id = NULL

scc = AddMetaData(scc,metadata = scibet)

source("~/genecorrect.R", echo=TRUE)
scc1 = gene.correct(scc)

saveRDS(scc1,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/SCC.rds')

#ccRCC---------------
post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/post.T.rds")
ccrcc = post.T
ccrcc$sample.timepoint = 'After'
ccrcc$treatment.type = if_else(ccrcc$donor_id %in% c('P915'),'PD1_CTLA4','PD1')
ccrcc$treatment.efficacy = if_else(ccrcc$ICB_Response =='ICB_PR','R','NR')
ccrcc$CancerType = 'ccRCC'
ccrcc$Study = 'ccRCC_Bi'
ccrcc$Platform = '10X'
ccrcc$sample.Tissue = 'Tumor'
ccrcc$DateType = 'Count'
ccrcc$sample.ID = paste0(ccrcc$Study,ccrcc$sample.timepoint,ccrcc$donor_id)

ccrcc@meta.data[,c(1,4:23)] = NULL

post.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/post.T.scibet.rds")
scibet = post.T.scibet
scibet$cell.id = NULL

ccrcc = AddMetaData(ccrcc,metadata = scibet)

source("~/genecorrect.R", echo=TRUE)
rcc = gene.correct(ccrcc)

saveRDS(rcc,file = "~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/ccRCC.rds")


# moshe ------------
pre.ctla <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_CTLA4_Moshe/pre.T.rds")
pre.dual <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/pre.T.rds")
post.dual <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/post.T.rds")
pre.pd1 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/pre.T.rds")
post.pd1 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/post.T.rds")

melan = merge(pre.ctla,c(pre.dual,post.dual,pre.pd1,post.pd1))

melan$sample.timepoint = if_else(melan$timepoint == 'Pre','Before','After')
melan$treatment.type = if_else(melan$therapy == 'anti-PD1','PD1',
                               if_else(melan$therapy == 'anti-CTLA4','CTLA4','PD1_CTLA4'))
melan$treatment.efficacy = if_else(melan$Response == 'Responder','R','NR')
melan$CancerType = 'Melanoma'
melan$Study = 'Melanoma_Feldman'
melan$Platform = 'SmartV2'
melan$sample.Tissue = 'Tumor'
melan$DateType = 'TPM'
melan$sample.ID = paste0(melan$Study,melan$sample.timepoint,melan$PtID)

melan@meta.data[,c(1,4:14)] = NULL

ctla4pre <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_CTLA4_Moshe/pre.T.scibet.rds")
dualpre <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/pre.T.scibet.rds")
dualpost <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/post.T.scibet.rds")
pd1pre <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/pre.T.scibet.rds")
pd1post <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/post.T.scibet.rds")

scibet = rbind(ctla4pre,dualpre,dualpost,pd1pre,pd1post)

scibet$cell.id = NULL

melan = AddMetaData(melan,scibet)

source("~/genecorrect.R", echo=TRUE)
melanm = gene.correct(melan)

saveRDS(melanm,file = "~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/Melanoma_Moshe.rds")


# livnat ------------------------------------------------------------------

post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/post.T.rds")
naive.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/naive.T.rds")

livnat = merge(post.T,naive.T)
livnat$sample.timepoint = if_else(livnat$treatment.group == 'post.treatment','After','NonTreat')
livnat$treatment.type = if_else(livnat$samples %in% c('Mel04.3','Mel58','Mel60','Mel88','Mel98'),'CTLA4',
                                if_else(livnat$samples %in% c('Mel121.1','Mel74'),'PD1',
                                        if_else(livnat$samples %in% c('Mel102','Mel106','Mel110','Mel126','Mel194','Mel72','Mel75',
                                                                      'Mel75.1','Mel78','Mel94'),'PD1_CTLA4','NonTreat')))
livnat$treatment.efficacy = 'NoInfo'
livnat$CancerType = 'Melanoma'
livnat$Study = 'Melanoma_Arnon'
livnat$Platform = 'SmartV2'
livnat$sample.Tissue = 'Tumor'
livnat$DateType = 'Count'
livnat$sample.ID = paste0(livnat$Study,livnat$sample.timepoint,livnat$samples)

livnat@meta.data[,c(1,4:10)] = NULL


post.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/post.T.scibet.rds")
naive.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/naive.T.scibet.rds")

scibet = rbind(post.T.scibet,naive.T.scibet)
scibet$cell.id = NULL

livnat = AddMetaData(livnat,scibet)

source("~/genecorrect.R", echo=TRUE)
livnat1 = gene.correct(livnat)
# do not have response info
saveRDS(livnat1,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/Melanoma_Livnat.rds')

# nsclc_liu---------

post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/post.T.rds")
pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/pre.T.rds")

liu1 = merge(pre.T,post.T)
liu$sample.timepoint = if_else(liu$timepoint == 'Pre','Before','After')
liu$treatment.type = if_else(liu$patient %in% c('P36','P37','P38','P1','P10','P13','P19','P29','P30','P33','P35'),'PD1','NonTreat')
liu$CancerType = 'NSCLC'
liu$Study = 'NSCLC_Liu'
liu$Platform = '10X'
liu$sample.Tissue = 'Tumor'
liu$DateType = 'Count'
liu$treatment.efficacy = if_else(liu$patient %in% c('P36','P37','P38'),'NR',
                                 if_else(liu$patient %in% c('P1','P10','P13','P19','P29','P30','P33','P35'),'R','NonTreat'))
liu$sample.ID = paste0(liu$Study,liu$sample.timepoint,liu$sample)

liu@meta.data[,c(1,4:9)] = NULL

pre.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/pre.T.scibet.rds")
post.T.scibet <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/post.T.scibet.rds")
scibet = rbind(pre.T.scibet,post.T.scibet)
scibet$cell.id = NULL

liu = AddMetaData(liu,scibet)

source("~/genecorrect.R", echo=TRUE)
liu = gene.correct(liu)

saveRDS(liu,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/NSCLC_Liu.rds')


# NSCLC_Caishu ------------------------------------------------------------
list.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Caushi/list.T.rds")

name.all = character()

for(i in 1:46){
  names = list.T[[i]]$sample[1]
  name.all[i] = names
}

caishu = merge(x = list.T[[1]],y = list.T[2:46],add.cell.ids = name.all)

caishu$sample.timepoint = 'After'
caishu$treatment.type = 'PD1'
caishu$CancerType = 'NSCLC'
caishu$Study = 'NSCLC_Caushi'
caishu$Platform = '10X'
caishu$sample.Tissue = 'Tumor'
caishu$DateType = 'Count'
caishu$sample.ID = paste0(caishu$Study,caishu$sample.timepoint,caishu$sample)

# MPR as R, non-MPR as NR
caishu$treatment.efficacy = if_else(caishu$sample %in% 
                                     c( "MD01-005_tumor_3","MD01-005_tumor_4","MD01-005_tumor_5", "MD01-005_tumor_6",    
                                        "MD01-005_tumor_7","MD01-005_tumor_8","MD01-005_tumor_9","MD01-010_tumor_1",
                                        "MD043-003_tumor_1","MD043-003_tumor_2","MD043-003_tumor_3","MD043-003_tumor_4",   
                                        "MD043-003_tumor_5","MD043-008_tumor_1","NY016-022_tumor_1","NY016-022_tumor_2",   
                                        "NY016-022_tumor_3","NY016-022_tumor_4","NY016-025_tumor_1","NY016-025_tumor_2",
                                        "NY016-025_tumor_3","NY016-025_tumor_4"),'R','NR')


caishu@meta.data[,c(1,4)] = NULL

scibet.list <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Caushi/scibet.list.rds")

data1 = scibet.list[[1]]
data1 = as.data.frame(data1)
rownames(data1) = paste0(names(scibet.list)[1],"_",data1$cell.id)

for (i in 2:length(scibet.list)) {
  data = scibet.list[[i]]
  data = as.data.frame(data)
  rownames(data) = paste0(names(scibet.list)[i],"_",data$cell.id)
  data1 = rbind(data1,data)
}

scibet = data1
scibet$cell.id = NULL

caishu = AddMetaData(caishu,scibet)

source("~/genecorrect.R", echo=TRUE)
caushi = gene.correct(caishu)

saveRDS(caushi,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Caushi/Caushi.rds')


# OSCC

new.list <- readRDS("~/PanCancerICI/Data/scDATA/OSCC/new.list.rds")
oscc = merge(new.list[[1]],y = new.list[2:25])
oscc = subset(oscc,CellType_ID == 'T cell')

oscc@meta.data[,c(1,4,9:11)] = NULL


oscc$sample.timepoint = if_else(oscc$Stage == 'Pre-Tx','Before','After')
oscc$treatment.type = if_else(oscc$Cohort == 'Mono','PD1','PD1_CTLA4')
oscc$CancerType = 'HNSCC'
oscc$Study = 'HNSCC_Luoma'
oscc$Platform = '10X'
oscc$sample.Tissue = 'Tumor'
oscc$DateType = 'Count'
oscc$sample.ID = paste0(oscc$Study,oscc$sample.timepoint,oscc$Patient_ID)
oscc$treatment.efficacy = if_else(oscc$Path_response == 'Low','NR','R')
oscc$scibet.celltype = 'notuse'
oscc@meta.data[,c(3:6)] = NULL
oscc1 = gene.correct(obj = oscc)

oscc1[["percent.mt"]] <- PercentageFeatureSet(oscc1, pattern = "^MT-")
oscc1$orig.ident = 'HNSCC'


# merge all -----------------------
BC_Bassez <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/BC_Bassez.rds")
BC_Zhang <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/BC_Zhang.rds")
BCC <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/BCC.rds")
ccRCC <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/ccRCC.rds")
Melanoma_Moshe <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/Melanoma_Moshe.rds")
Melanoma_Livnat <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/Melanoma_Livnat.rds")
Caushi <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Caushi/Caushi.rds")
NSCLC_Liu <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/NSCLC_Liu.rds")
SCC <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/SCC.rds")
HNSCC = readRDS("~/HNSCC.T.rds")

# Caushi data have multiple sample for one pt
Caushi = subset(Caushi,sample.ID %in% c('NSCLC_CaushiAfterMD01-004_tumor_1','NSCLC_CaushiAfterMD01-005_tumor_7','NSCLC_CaushiAfterMD01-010_tumor_1',
                                        'NSCLC_CaushiAfterMD01-019_tumor_1','NSCLC_CaushiAfterMD01-024_tumor_1','NSCLC_CaushiAfterMD043-003_tumor_3',
                                        'NSCLC_CaushiAfterMD043-006_tumor_1','NSCLC_CaushiAfterMD043-008_tumor_1','NSCLC_CaushiAfterMD043-011_tumor_3',
                                        'NSCLC_CaushiAfterNY016-007_tumor_1','NSCLC_CaushiAfterNY016-014_tumor_1','NSCLC_CaushiAfterNY016-015_tumor_2',
                                        'NSCLC_CaushiAfterNY016-021_tumor_1','NSCLC_CaushiAfterNY016-022_tumor_1','NSCLC_CaushiAfterNY016-025_tumor_2'))


# Liu data have multiple sample for one pt
NSCLC_Liu = subset(NSCLC_Liu,sample.ID %in% c(names(table(NSCLC_Liu$sample.ID))[c(-1,-3,-6)]) )

all = merge(BC_Bassez,c(BC_Zhang,BCC,Caushi,ccRCC,Melanoma_Livnat,Melanoma_Moshe,NSCLC_Liu,SCC,OSCC),
            add.cell.id = c('BC_Bassez','BC_Zhang','BCC','Caushi','ccRCC','Melanoma_Livnat','Melanoma_Moshe','NSCLC_Liu','SCC','HNSCC'))

table(all$Study,all$sample.timepoint)
table(all$Study,all$treatment.type)
table(all$Study,all$treatment.efficacy)
table(all$Study,all$Platform)
table(all$Study,all$DateType)
table(all$Study,all$sample.Tissue)
table(all$Study,all$sample.timepoint)

all = subset(all,sample.timepoint != 'NonTreat')# filter 1438 Melanoma_Arnon
all = subset(all,treatment.type != 'NonTreat')# filter 56103 liu
all = subset(all,treatment.efficacy != 'NoInfo')# filter 1883 Melanoma_Arnon
all = subset(all,DateType != 'TPM') # filter 8436 Melanoma_Feldman

all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^MT-")
all = subset(all,percent.mt<20) # 
VlnPlot(all,features = 'percent.mt',split.by = 'Study',pt.size = 0)

saveRDS(all,file = 'all.T.filter.mt.rds')



# geneblacklist-----------
source("~/Final/geneblacklist.R", echo=TRUE)
rpgene = grep('^RP[LS]',rownames(HNSCC.T),value = T)
mtgene = 'MALAT1'
final.gene = setdiff(rownames(HNSCC.T),y = c(immgene,rpgene,mtgene))
HNSCC.T1 = HNSCC.T[final.gene,]
saveRDS(all1,file = 'all.T.filter.mt.geneblacklist.rds')

# remove HSP gene
all.T.filter.mt.geneblacklist <- readRDS("~/all.T.filter.mt.geneblacklist.rds")
HSPgene = grep('^HSP',rownames(all.T.filter.mt.geneblacklist),value = T)
final.gene = setdiff(rownames(all.T.filter.mt.geneblacklist),HSPgene)
all2 = all.T.filter.mt.geneblacklist[final.gene,]
saveRDS(all2,file = 'all.T.filter.mt.geneblacklist.hsp.rds')

# filter <3 gene, <200 cell
count = all2[['RNA']]@counts
seu = CreateSeuratObject(counts = count,min.cells = 3,min.features = 100,meta.data = all.T.filter.mt.geneblacklist@meta.data)
saveRDS(seu,file = 'all.T.filter.mt.geneblacklist.hsp.minC3F100.rds')

# kegg--------
#kegg = GSEABase::getGmt('~/PanCancerICI/CD8TMetabolism/Data/KEGG_metabolism_nc.gmt')

#kegg.list = list()
#for (i in 1:length(kegg)) {
#  kegg.list[[kegg[[i]]@setName]] =  kegg[[i]]@geneIds
#}

#kegg.gene = kegg.list[[1]]

#for (i in 2:length(kegg.list)) {
#  gene = kegg.list[[i]]
#  kegg.gene = append(kegg.gene,gene)
#}

#kegg.gene = kegg.gene[!duplicated(kegg.gene)]

kegg.gene = readRDS('~/PanCancerICI/CD8TMetabolism/Data/metabolic.gene.name.rds')
all.kegg = seu[kegg.gene,]
saveRDS(all.kegg,file = 'all.T.filter.mt.geneblacklist.hsp.minC3F100.kegg.rds')


