## bulkDATA standardized

# def: expr = TPM/FPKN, count = rawcount
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
# 1 Braun_2020_ccRCC_PD1 -----------------------------------------------------------------
# ref: Interplay of somatic alterations and immune infiltration modulates response to PD-1 blockade in advanced clear cell renal cell carcinoma

normed <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/ccRCC/normed.rds")
clinical <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/ccRCC/clinical.rds")

sample_n = colnames(normed)
all_n = clinical$RNA_ID
all_n = all_n[all_n!='NA']
table(sample_n %in% all_n)

clinical_subset = dplyr::filter(clinical,RNA_ID %in% all_n)
clinical_subset = dplyr::filter(clinical_subset, Arm == 'NIVOLUMAB')

normed = normed[!duplicated(normed$gene_name),]

normed_subset = normed[,colnames(normed) %in% clinical_subset$RNA_ID]

Braun_2020_ccRCC_PD1 = list(clinical = clinical_subset, expr = normed_subset)

saveRDS(Braun_2020_ccRCC_PD1,'Braun_2020_ccRCC_PD1.rds')

expr = BulkICIdata$Braun_2020_ccRCC_PD1$expr
expr$symbol = normed$gene_name

expr = expr[!duplicated(expr$symbol),]
expr = as.data.frame(expr)
rownames(expr) = expr$symbol
expr$symbol = NULL

BulkICIdata$Braun_2020_ccRCC_PD1$expr = expr

# 2 Cho_2020_NSCLC_PD1 ------------------------------------------------------------
# ref: Genome-wide identification of differentially methylated promoters and enhancers associated with response to anti-PD-1 therapy in non-small cell lung cancer

count <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/GSE126043_NSCLC/count.rds")
clinical <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/GSE126043_NSCLC/clinical.rds")

Cho_2020_NSCLC_PD1 = list(clinical = clinical, count = count)
saveRDS(Cho_2020_NSCLC_PD1,'Cho_2020_NSCLC_PD1.rds')

count = BulkICIdata$Cho_2020_NSCLC_PD1$count
table(duplicated(count$gene))
count = count[!duplicated(count$gene),]
rownames(count) = count$gene
count$gene = NULL

BulkICIdata$Cho_2020_NSCLC_PD1$count = count

# 3 Jung_2019_NSCLC_PD1_PDL1 ----------------------------------------------------------------
# ref: DNA methylation loss promotes immune evasion of tumours with high mutation and copy number load
# note: DCB = benefit NCB = nonbenefit

tpm <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/GSE135222_NSCLC/tpm.rds")
PFS <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/GSE135222_NSCLC/PFS.rds")

Jung_2019_NSCLC_PD1_PDL1 = list(clinical = PFS, expr = tpm)
saveRDS(Jung_2019_NSCLC_PD1_PDL1,'Jung_2019_NSCLC_PD1_PDL1.rds')

tpm = BulkICIdata$Jung_2019_NSCLC_PD1_PDL1$expr
tpm$ENSEMBL<-substring(tpm$gene_id, 1, 15)
geneid<-bitr(tpm$ENSEMBL, fromType='ENSEMBL', toType='SYMBOL', OrgDb='org.Hs.eg.db', drop = TRUE)
tpm = left_join(tpm,geneid,by = 'ENSEMBL')
tpm = tpm[!is.na(tpm$SYMBOL),]
tpm <- tpm[!duplicated(tpm$SYMBOL),]
rownames(tpm) = tpm$SYMBOL
tpm$SYMBOL = NULL
tpm$gene_id = NULL
tpm$ENSEMBL = NULL
BulkICIdata$Jung_2019_NSCLC_PD1_PDL1$expr = tpm

# 4 Hsu_2021_HCC_PD1_PDL1 -------------------------------------------------
# ref: Exploring Markers of Exhausted CD8 T Cells to Predict Response to Immune Checkpoint Inhibitor Therapy for Hepatocellular Carcinoma

tpm <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/GSE141016_HCC/tpm.rds")
clinical <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/GSE141016_HCC/clinical.rds")

Hsu_2021_HCC_PD1_PDL1 = list(clinical = clinical, expr = tpm)
saveRDS(Hsu_2021_HCC_PD1_PDL1,'Hsu_2021_HCC_PD1_PDL1.rds')

tpm = BulkICIdata$Hsu_2021_HCC_PD1_PDL1$expr
table(duplicated(tpm$ID_REF))
rownames(tpm) = tpm$ID_REF
tpm$ID_REF = NULL
BulkICIdata$Hsu_2021_HCC_PD1_PDL1$expr = tpm

# 5 Hugo_2016_Melanoma_PD1 ------------------------------------------------
# ref: Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma

tpm <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/GSE78220_melanoma/tpm.rds")
respon <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/GSE78220_melanoma/respon.rds")
os <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/GSE78220_melanoma/os.rds")

Hugo_2016_Melanoma_PD1 = list(clinical = os, expr = tpm)
saveRDS(Hugo_2016_Melanoma_PD1,'Hugo_2016_Melanoma_PD1.rds')

tpm = BulkICIdata$Hugo_2016_Melanoma_PD1$expr
table(duplicated(tpm$Gene))
tpm = as.data.frame(tpm)
rownames(tpm) = tpm$Gene
tpm$Gene = NULL
BulkICIdata$Hugo_2016_Melanoma_PD1$expr = tpm

# 6 Riaz_2017_Melanoma_PD1 ------------------------------------------------
# ref: Tumor and Microenvironment Evolution during Immunotherapy with Nivolumab

count <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/GSE91061_melanoma/count.rds")
clinical <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/GSE91061_melanoma/clinical.rds")

name = vector()
for (i in 2:(length(colnames(count))-1)) {
  if(strsplit(colnames(count),'_')[[i]][2] == 'Pre'){
    name[i] = 'pre'
  }else if(strsplit(colnames(count),'_')[[i]][2] == 'On'){
    name[i] = 'on'
  }else {
    name[i] = colnames(count)[i]
  }
}
name[1] = 'ENTRIEZID'
name[111] = 'SYMBOL'
pre.expr = count[,name =='pre']
pre.expr$SYMBOL = count$SYMBOL
on.expr = count[,name =='on']
on.expr$SYMBOL = count$SYMBOL

count.list = list(pre.expr,on.expr)

Riaz_2017_Melanoma_PD1 = list(clinical = clinical, count = count.list)
saveRDS(Riaz_2017_Melanoma_PD1,file = 'Riaz_2017_Melanoma_PD1.rds')

count = BulkICIdata$Riaz_2017_Melanoma_PD1$count
table(duplicated(count[[1]]$SYMBOL))
count[[1]] = count[[1]][!duplicated(count[[1]]$SYMBOL),]
count[[1]] = count[[1]][!is.na(count[[1]]$SYMBOL),]
rownames(count[[1]]) = count[[1]]$SYMBOL
count[[1]]$SYMBOL = NULL

table(duplicated(count[[2]]$SYMBOL))
count[[2]] = count[[2]][!duplicated(count[[2]]$SYMBOL),]
count[[2]] = count[[2]][!is.na(count[[2]]$SYMBOL),]
rownames(count[[2]]) = count[[2]]$SYMBOL
count[[2]]$SYMBOL = NULL

names(count) = c('Pre','On')
BulkICIdata$Riaz_2017_Melanoma_PD1$count = count
count1 = BulkICIdata$Riaz_2017_Melanoma_PD1$count[[1]]
count1 = count1[!duplicated(count1$SYMBOL),]
count2 = BulkICIdata$Riaz_2017_Melanoma_PD1$count[[2]]

BulkICIdata$Riaz_2017_Melanoma_PD1$count[[1]] = BulkICIdata$Riaz_2017_Melanoma_PD1$count[[1]][!is.na(BulkICIdata$Riaz_2017_Melanoma_PD1$count[[1]]$SYMBOL),]
rownames(BulkICIdata$Riaz_2017_Melanoma_PD1$count[[1]]) = BulkICIdata$Riaz_2017_Melanoma_PD1$count[[1]]$SYMBOL
BulkICIdata$Riaz_2017_Melanoma_PD1$count[[1]]$SYMBOL = NULL
BulkICIdata$Riaz_2017_Melanoma_PD1$count[[2]] = BulkICIdata$Riaz_2017_Melanoma_PD1$count[[2]][!duplicated(BulkICIdata$Riaz_2017_Melanoma_PD1$count[[2]]$SYMBOL),]
BulkICIdata$Riaz_2017_Melanoma_PD1$count[[2]] = BulkICIdata$Riaz_2017_Melanoma_PD1$count[[2]][!is.na(BulkICIdata$Riaz_2017_Melanoma_PD1$count[[2]]$SYMBOL),]
rownames(BulkICIdata$Riaz_2017_Melanoma_PD1$count[[2]]) = BulkICIdata$Riaz_2017_Melanoma_PD1$count[[2]]$SYMBOL
BulkICIdata$Riaz_2017_Melanoma_PD1$count[[2]]$SYMBOL = NULL


# 7 Mariathasan_2018_UC_PDL1 ---------------------------------------------
# ref: TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells
# note: IMvigor210

clinical <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/IMvigor210_bladder/clinical.rds")
count <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/IMvigor210_bladder/count.rds")

Mariathasan_2018_UC_PDL1 = list(clinical = clinical,count = count)
saveRDS(Mariathasan_2018_UC_PDL1,'Mariathasan_2018_UC_PDL1.rds')


# 8 Liu_2019_Melanoma_PD1 -------------------------------------------------
# ref: Integrative molecular and clinical modeling of clinical outcomes to PD1 blockade in patients with metastatic melanoma

clinical <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/meta_melanoma/clinical.rds")
tpm <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/meta_melanoma/tpm.rds")

Liu_2019_melanoma_PD1 = list(clinical = clinical, expr = tpm)
saveRDS(Liu_2019_melanoma_PD1,'Liu_2019_Melanoma_PD1.rds')
expr = BulkICIdata$Liu_2019_Melanoma_PD1$expr
rownames(expr) = gsub('\\.','-', rownames(expr))
BulkICIdata$Liu_2019_Melanoma_PD1$expr = expr

# 9 Gide_2019_Melanoma_PD1_CTLA4 ----------------------------------------------------
# ref: Distinct Immune Cell Populations Define Response to Anti-PD-1 Monotherapy and Anti-PD-1/Anti-CTLA-4 Combined Therapy

count <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/PRJEB23709_melanoma/count.rds")
clinical_pd1 <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/PRJEB23709_melanoma/clinical_pd1.rds")
clinical_ipipd1 <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/PRJEB23709_melanoma/clinical_ipipd1.rds")

Gide_2019_Melanoma_PD1_CTLA4 = list(clinical = list(Comb = clinical_ipipd1, Mono = clinical_pd1),count = count)
saveRDS(Gide_2019_Melanoma_PD1_CTLA4,'Gide_2019_Melanoma_PD1_CTLA4.rds')

count = BulkICIdata$Gide_2019_Melanoma_PD1_CTLA4$count
geneid<-bitr(count$ID, fromType='ENSEMBL', toType='SYMBOL', OrgDb='org.Hs.eg.db', drop = TRUE)
colnames(count)[1] ='ENSEMBL'
count = left_join(count,geneid,by = 'ENSEMBL')
count = count[!is.na(count$SYMBOL),]
count <- count[!duplicated(count$SYMBOL),]
rownames(count) = count$SYMBOL
count$SYMBOL = NULL
count$ENSEMBL = NULL

BulkICIdata$Gide_2019_Melanoma_PD1_CTLA4$count =count

# 10 Zhao_2019_GBM_PD1 ----------------------------------------------------
# ref: Immune and genomic correlates of response to anti-PD-1 immunotherapy in glioblastoma

tpm <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/PRJNA482620_glioma/tpm.rds")
count <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/PRJNA482620_glioma/count.rds")
clinical <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/PRJNA482620_glioma/clinical.rds")
# adopted from ZZ paper IOS

clinical = ios$Zhao_GBM_pre_aPD1[,(1:7)]
clinical_meta = ls_meta$Zhao_GBM_pre_aPD1

count = BulkICIdata$Zhao_2019_GBM_PD1$count
count$ENSEMBL<-substring(rownames(count), 1, 15)
geneid<-bitr(count$ENSEMBL, fromType='ENSEMBL', toType='SYMBOL', OrgDb='org.Hs.eg.db', drop = TRUE)
count = left_join(count,geneid,by = 'ENSEMBL')
count = count[!is.na(count$SYMBOL),]
count <- count[!duplicated(count$SYMBOL),]
rownames(count) = count$SYMBOL
count$SYMBOL = NULL
count$ENSEMBL = NULL

tpm = BulkICIdata$Zhao_2019_GBM_PD1$expr
tpm$ENSEMBL<-substring(rownames(tpm), 1, 15)
geneid<-bitr(tpm$ENSEMBL, fromType='ENSEMBL', toType='SYMBOL', OrgDb='org.Hs.eg.db', drop = TRUE)
tpm = left_join(tpm,geneid,by = 'ENSEMBL')
tpm = tpm[!is.na(tpm$SYMBOL),]
tpm <- tpm[!duplicated(tpm$SYMBOL),]
rownames(tpm) = tpm$SYMBOL
tpm$SYMBOL = NULL
tpm$ENSEMBL = NULL


Zhao_2019_GBM_PD1 = list(clinical = clinical, count = count,expr = tpm)
saveRDS(Zhao_2019_GBM_PD1,'Zhao_2019_GBM_PD1.rds')

BulkICIdata$Zhao_2019_GBM_PD1$count = count
BulkICIdata$Zhao_2019_GBM_PD1$expr = tpm

# 11 Snyder_2017_UC_PD1 ---------------------------------------------------
# ref: Contribution of systemic and somatic factors to clinical response and resistance to PD‑L1 blockade in urothelial cancer: an exploratory multi‑omic analysis.

count <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/urothelial/count.rds")
clinical <- readRDS("~/Project/PanCancerICI/Data/bulkDATA/urothelial/clinical.rds")
data_pdl1_all <- read_csv("Project/PanCancerICI/Data/bulkDATA/urothelial/multi-omic-urothelial-anti-pdl1-master/data_pdl1_all.csv")

meta = data_pdl1_all[,c(3:11,13)]
count_id = colnames(count)[-1]
count_id = substr(count_id,3,100)
for (i in 1:25) {
  if(nchar(count_id[i])==2){
    count_id[i] = paste0('00',count_id[i])
  }else if(nchar(count_id[i])==3){
    count_id[i] = paste0('0',count_id[i])
  }else{count_id[i] = count_id[i]}
}
meta = dplyr::filter(meta,patient_id %in% count_id)

Snyder_2017_UC_PD1 = list(clinical = meta, count = count)
saveRDS(Snyder_2017_UC_PD1,'Snyder_2017_UC_PD1.rds')

count = BulkICIdata$Snyder_2017_UC_PD1$count
table(duplicated(count$ptgene_name))
count = count[!duplicated(count$ptgene_name),]
count = count[!is.na(count$ptgene_name),]
rownames(count) = count$ptgene_name
count$ptgene_name = NULL
BulkICIdata$Snyder_2017_UC_PD1$count = count


# 12 Kim_2018_GC_PD1 ------------------------------------------------------
# ref: Comprehensive molecular characterization of clinical responses to PD‑1 inhibition in metastatic gastric cancer
# note: original data unavailable, adopted from zz paper ios

clinical = ios$Kim_GC_pre_aPD1[,c(1,6)] # 1 = R, 0 = NR
expr = ios$Kim_GC_pre_aPD1[,c(1,8:34012)]

Kim_2018_GC_PD1 = list(clinical = clinical, expr = expr)
saveRDS(Kim_2018_GC_PD1,'Kim_2018_GC_PD1.rds')

tpm = BulkICIdata$Kim_2018_GC_PD1$expr
tpm = t(tpm)
tpm = as.data.frame(tpm)
tpm= tpm[-1,]
BulkICIdata$Kim_2018_GC_PD1$expr = tpm

# 13 Van_2015_Melanoma_CTLA4 ------------------------------------------------------------------
# ref: Genomic correlates of response to CTLA-4 blockade in metastatic melanoma
# note: original data unavailable, adopted from zz paper ios

clinical = ios$Van_SKCM_pre_aPD1[,c(1:6)] # 1 = R, 0 = NR
meta = ls_meta$Van_SKCM_pre_aPD1
clinical = dplyr::left_join(clinical,meta, by='ID')
expr = ios$Kim_GC_pre_aPD1[,c(1,8:25664)]

Van_2015_Melanoma_CTLA4 = list(clinical = clinical, expr = expr)
saveRDS(Van_2015_Melanoma_CTLA4,'Van_2015_Melanoma_CTLA4.rds')

tpm = BulkICIdata$Van_2015_Melanoma_CTLA4$expr
tpm = t(tpm)
tpm = as.data.frame(tpm)
tpm= tpm[-1,]
BulkICIdata$Van_2015_Melanoma_CTLA4$expr = tpm

# Merge -------------------------------------------------------------------

Van_2015_Melanoma_CTLA4 <- readRDS("~/Van_2015_Melanoma_CTLA4.rds")
Kim_2018_GC_PD1 <- readRDS("~/Kim_2018_GC_PD1.rds")
Snyder_2017_UC_PD1 <- readRDS("~/Snyder_2017_UC_PD1.rds")
Zhao_2019_GBM_PD1 <- readRDS("~/Zhao_2019_GBM_PD1.rds")
Gide_2019_Melanoma_PD1_CTLA4 <- readRDS("~/Gide_2019_Melanoma_PD1_CTLA4.rds")
Liu_2019_Melanoma_PD1 <- readRDS("~/Liu_2019_Melanoma_PD1.rds")
Mariathasan_2018_UC_PDL1 <- readRDS("~/Mariathasan_2018_UC_PDL1.rds")
Riaz_2017_Melanoma_PD1 <- readRDS("~/Riaz_2017_Melanoma_PD1.rds")
Hugo_2016_Melanoma_PD1 <- readRDS("~/Hugo_2016_Melanoma_PD1.rds")
Hsu_2021_HCC_PD1_PDL1 <- readRDS("~/Hsu_2021_HCC_PD1_PDL1.rds")
Jung_2019_NSCLC_PD1_PDL1 <- readRDS("~/Jung_2019_NSCLC_PD1_PDL1.rds")
Cho_2020_NSCLC_PD1 <- readRDS("~/Cho_2020_NSCLC_PD1.rds")
Braun_2020_ccRCC_PD1 <- readRDS("~/Braun_2020_ccRCC_PD1.rds")

BulkICIdata = list(Van_2015_Melanoma_CTLA4 = Van_2015_Melanoma_CTLA4,
                   Hugo_2016_Melanoma_PD1 = Hugo_2016_Melanoma_PD1,
                   Riaz_2017_Melanoma_PD1 = Riaz_2017_Melanoma_PD1,
                   Snyder_2017_UC_PD1 = Snyder_2017_UC_PD1,
                   Kim_2018_GC_PD1 = Kim_2018_GC_PD1,
                   Mariathasan_2018_UC_PDL1 = Mariathasan_2018_UC_PDL1,
                   Gide_2019_Melanoma_PD1_CTLA4 = Gide_2019_Melanoma_PD1_CTLA4,
                   Liu_2019_Melanoma_PD1 = Liu_2019_Melanoma_PD1,
                   Jung_2019_NSCLC_PD1_PDL1 = Jung_2019_NSCLC_PD1_PDL1,
                   Zhao_2019_GBM_PD1 = Zhao_2019_GBM_PD1,
                   Cho_2020_NSCLC_PD1 = Cho_2020_NSCLC_PD1,
                   Braun_2020_ccRCC_PD1 = Braun_2020_ccRCC_PD1,
                   Hsu_2021_HCC_PD1_PDL1 = Hsu_2021_HCC_PD1_PDL1)
saveRDS(BulkICIdata,'BulkICIdata.rds')


# load data ---------------------------------------------------------------

BulkICIdata <- readRDS("~/Project/Signatures/Bulk_ICI_data/BulkICIdata.rds")

status = vector()
for (i in 1:length(BulkICIdata)) {
  testdata = BulkICIdata[[i]]
  obj = names(testdata)
  if (('count' %in% obj) & ('expr' %in% obj)) {
    status[i] = c('Finished')
  } else if (('count' %in% obj) & !('expr' %in% obj)) {
    status[i] = c('Count only')  # calculate TPM use IOBR
  } else if (!('count' %in% obj) & ('expr' %in% obj)) {
    status[i] = c('Expr only')
    testdata$count = 'Unavailable'
  }
}

# status 3, 4, 6, 7, 11 count---------------
"""
data = BulkICIdata[[3]]
names(data$count) = c('Pre','On')

countpre = data$count[[1]]
tpm_pre = count2tpm(
    countMat = countpre,
    idType = 'SYMBOL',
    org = 'hsa',
    gene_symbol = rownames(countpre)
  )
counton = data$count[[2]]
tpm_on = count2tpm(
  countMat = counton,
  idType = 'SYMBOL',
  org = 'hsa',
  gene_symbol = rownames(counton)
)

data$expr = list(Pre = tpm_pre, On = tpm_on)
BulkICIdata[[3]] = data

data = BulkICIdata[[4]]
count = data$count
tpm = count2tpm(
  countMat = count,
  idType = 'SYMBOL',
  org = 'hsa',
  gene_symbol = rownames(count)
)

data$expr = tpm
BulkICIdata[[4]] = data


data = BulkICIdata[[6]]
count = data$count
tpm = count2tpm(
  countMat = count,
  idType = 'SYMBOL',
  org = 'hsa',
  gene_symbol = rownames(count)
)
data$expr = tpm
BulkICIdata[[6]] = data

data = BulkICIdata[[7]]
count = data$count
tpm = count2tpm(
  countMat = count,
  idType = 'SYMBOL',
  org = 'hsa',
  gene_symbol = rownames(count)
)
data$expr = tpm
BulkICIdata[[7]] = data

data = BulkICIdata[[11]]
count = data$count
tpm = count2tpm(
  countMat = count,
  idType = 'SYMBOL',
  org = 'hsa',
  gene_symbol = rownames(count)
)
data$expr = tpm
BulkICIdata[[11]] = data
saveRDS(BulkICIdata,'BulkICIdata.rds')
"""


# correct gene ------------------------------------------------------------
source("~/Project/PaperCD8/code/Bulk/geneCorrect.bulk.R")
BulkICIdata <- readRDS("~/Project/Signatures/Bulk_ICI_data/BulkICIdata.rds")

for (i in 1:length(BulkICIdata)) {
  data = BulkICIdata[[i]]$expr
  if (i == 3) {
    data1 = data$Pre
    new1 = gene.correct(data1)
    BulkICIdata[[i]]$expr$Pre = new1
    data2 = data$On
    new2 = gene.correct(data2)
    BulkICIdata[[i]]$expr$On = new2
  } else {
    new = gene.correct(data)
    BulkICIdata[[i]]$expr = new
  }
}

for (i in c(3, 4, 6, 7, 11)) {
  data = BulkICIdata[[i]]$count
  if (i == 3) {
    data1 = data$Pre
    new1 = gene.correct(data1)
    BulkICIdata[[i]]$count$Pre = new1
    data2 = data$On
    new2 = gene.correct(data2)
    BulkICIdata[[i]]$count$On = new2
  } else {
    new = gene.correct(data)
    BulkICIdata[[i]]$count = new
  }
}
