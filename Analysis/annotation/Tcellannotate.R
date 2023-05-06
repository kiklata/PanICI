run.singler = function(obj){
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(BiocParallel)
  library(SingleCellExperiment)
  
  ref = celldex::MonacoImmuneData()
  seu = obj
  seu = as.SingleCellExperiment(seu)
  seu.singler = SingleR(test = seu, ref = ref, labels = ref$label.main, fine.tune = T,
                        de.method="wilcox",de.n = 50, BPPARAM=MulticoreParam(16))
  return(seu.singler)
}

# Melanoma_Moshe---------
#PD1
pre <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/pre.rds")
post <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/post.rds")

#pre
seu.singler = run.singler(pre)

saveRDS(seu.singler,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/pre.singler.rds')
table(seu.singler$pruned.labels)

pre$singler = seu.singler$pruned.labels
seu.t = subset(pre,singler %in% c('CD4+ T cells','CD8+ T cells','T cells'))
saveRDS(seu.t,"~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/pre.T.rds")

#post
seu.singler = run.singler(post)

saveRDS(seu.singler,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/post.singler.rds')
table(seu.singler$pruned.labels)

post$singler = seu.singler$pruned.labels
seu.t = subset(post,singler %in% c('CD4+ T cells','CD8+ T cells','T cells'))
saveRDS(seu.t,"~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/post.T.rds")

#PD1+CTLA4
pre <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/pre.rds")
post <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/post.rds")

#pre
seu.singler = run.singler(pre)

saveRDS(seu.singler,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/pre.singler.rds')
table(seu.singler$pruned.labels)

pre$singler = seu.singler$pruned.labels
seu.t = subset(pre,singler %in% c('CD4+ T cells','CD8+ T cells','T cells'))
saveRDS(seu.t,"~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/pre.T.rds")

#post
seu.singler = run.singler(post)

saveRDS(seu.singler,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/post.singler.rds')
table(seu.singler$pruned.labels)

post$singler = seu.singler$pruned.labels
seu.t = subset(post,singler %in% c('CD4+ T cells','CD8+ T cells','T cells'))
saveRDS(seu.t,"~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/post.T.rds")

# CTLA4 
pre <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_CTLA4_Moshe/pre.rds")
seu.singler = run.singler(pre)

saveRDS(seu.singler,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_CTLA4_Moshe/pre.singler.rds')
table(seu.singler$pruned.labels)

pre$singler = seu.singler$pruned.labels
seu.t = subset(pre,singler %in% c('CD4+ T cells','CD8+ T cells','T cells'))
saveRDS(seu.t,"~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_CTLA4_Moshe/pre.T.rds")

#ccrcc--------------
post <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/post.rds")
seu.singler = run.singler(post)
saveRDS(seu.singler,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/post.singler.rds')
table(seu.singler$pruned.labels)

post$singler = seu.singler$pruned.labels
seu.t = subset(post,singler %in% c('CD4+ T cells','CD8+ T cells','T cells'))
saveRDS(seu.t,"~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/post.T.rds")


# run SciBet--------------------

run.scibet = function(obj,ref){
  
  library(scibet)
  library(Seurat)
  
  expr = obj[['RNA']]@counts
  expr = scuttle::calculateTPM(expr)
  expr = as.matrix(expr)
  expr = t(expr)
  prd <- SciBet(ref, expr)
  prd = as.data.frame(prd)
  prd$cellid = rownames(expr)
  
  colnames(prd) = c('scibet.celltype','cell.id')
  rownames(prd) = prd$cell.id
  
  return(prd)
}

SciBet.ref <- readRDS("~/software/scibet/count/SciBet.ref.rds")

CD8subset = c("CD8.c01.Tn.MAL",        "CD8.c02.Tm.IL7R",       "CD8.c03.Tm.RPS12",     
              "CD8.c04.Tm.CD52",       "CD8.c05.Tem.CXCR5",     "CD8.c06.Tem.GZMK",     
              "CD8.c07.Temra.CX3CR1",  "CD8.c08.Tk.TYROBP",     "CD8.c09.Tk.KIR2DL4",   
              "CD8.c10.Trm.ZNF683",    "CD8.c11.Tex.PDCD1",     "CD8.c12.Tex.CXCL13",   
              "CD8.c13.Tex.myl12a",    "CD8.c14.Tex.TCF7",      "CD8.c15.ISG.IFIT1",    
              "CD8.c16.MAIT.SLC4A10",  "CD8.c17.Tm.NME1")

#BC_PD1---------------
cohort1.pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.pre.T.rds")
cohort1.pre.T.scibet = run.scibet(obj = cohort1.pre.T,ref = SciBet.ref)
saveRDS(cohort1.pre.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.pre.T.scibet.rds')
cohort1.pre.T = AddMetaData(cohort1.pre.T,cohort1.pre.T.scibet)
cohort1.pre.T.cd8 = subset(cohort1.pre.T,scibet.celltype %in% CD8subset)
saveRDS(cohort1.pre.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.pre.T.cd8.rds')

cohort1.on.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.on.T.rds")
cohort1.on.T.scibet = run.scibet(obj = cohort1.on.T,ref = SciBet.ref)
saveRDS(cohort1.on.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.on.T.scibet.rds')
cohort1.on.T = AddMetaData(cohort1.on.T,cohort1.on.T.scibet)
cohort1.on.T.cd8 = subset(cohort1.on.T,scibet.celltype %in% CD8subset)
saveRDS(cohort1.on.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort1.on.T.cd8.rds')

cohort2.pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.pre.T.rds")
cohort2.pre.T.scibet = run.scibet(obj = cohort2.pre.T,ref = SciBet.ref)
saveRDS(cohort2.pre.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.pre.T.scibet.rds')
cohort2.pre.T = AddMetaData(cohort2.pre.T,cohort2.pre.T.scibet)
cohort2.pre.T.cd8 = subset(cohort2.pre.T,scibet.celltype %in% CD8subset)
saveRDS(cohort2.pre.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.pre.T.cd8.rds')

cohort2.on.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.on.T.rds")
cohort2.on.T.scibet = run.scibet(obj = cohort2.on.T,ref = SciBet.ref)
saveRDS(cohort2.on.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.on.T.scibet.rds')
cohort2.on.T = AddMetaData(cohort2.on.T,cohort2.on.T.scibet)
cohort2.on.T.cd8 = subset(cohort2.on.T,scibet.celltype %in% CD8subset)
saveRDS(cohort2.on.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PD1/cohort2.on.T.cd8.rds')

#BC_PDL1----------------

pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/pre.T.rds")
pre.T.scibet = run.scibet(obj = pre.T,ref = SciBet.ref)
saveRDS(pre.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/pre.T.scibet.rds')
pre.T = AddMetaData(pre.T,pre.T.scibet)
pre.T.cd8 = subset(pre.T,scibet.celltype %in% CD8subset)
saveRDS(pre.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/pre.T.cd8.rds')

post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/post.T.rds")
post.T.scibet = run.scibet(obj = post.T,ref = SciBet.ref)
saveRDS(post.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/post.T.scibet.rds')
post.T = AddMetaData(post.T,post.T.scibet)
post.T.cd8 = subset(post.T,scibet.celltype %in% CD8subset)
saveRDS(post.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BC_PDL1/post.T.cd8.rds')

#BCC_PD1--------------

pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/pre.T.rds")
pre.T.scibet = run.scibet(obj = pre.T,ref = SciBet.ref)
saveRDS(pre.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/pre.T.scibet.rds')
pre.T = AddMetaData(pre.T,pre.T.scibet)
pre.T.cd8 = subset(pre.T,scibet.celltype %in% CD8subset)
saveRDS(pre.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/pre.T.cd8.rds')

post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/post.T.rds")
post.T.scibet = run.scibet(obj = post.T,ref = SciBet.ref)
saveRDS(post.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/post.T.scibet.rds')
post.T = AddMetaData(post.T,post.T.scibet)
post.T.cd8 = subset(post.T,scibet.celltype %in% CD8subset)
saveRDS(post.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/BCC_PD1/post.T.cd8.rds')

#ccRCC---------------
#annotate

post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/post.T.rds")
post.T.scibet = run.scibet(obj = post.T,ref = SciBet.ref)
saveRDS(post.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/post.T.scibet.rds')
post.T = AddMetaData(post.T,post.T.scibet)
post.T.cd8 = subset(post.T,scibet.celltype %in% CD8subset)
saveRDS(post.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/ccRCC_Kevin/post.T.cd8.rds')

#Melanoma_Moshe_ctla4------------------
# TPM data only, thus run.scibet() here required no TPM calc
# !!

pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_CTLA4_Moshe/pre.T.rds")
pre.T.scibet = run.scibet(obj = pre.T,ref = SciBet.ref)
saveRDS(pre.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_CTLA4_Moshe/pre.T.scibet.rds')
pre.T = AddMetaData(pre.T,pre.T.scibet)
pre.T.cd8 = subset(pre.T,scibet.celltype %in% CD8subset)
saveRDS(pre.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_CTLA4_Moshe/pre.T.cd8.rds')

#pd1_ctla4-----------------
pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/pre.T.rds")
pre.T.scibet = run.scibet(obj = pre.T,ref = SciBet.ref)
saveRDS(pre.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/pre.T.scibet.rds')
pre.T = AddMetaData(pre.T,pre.T.scibet)
pre.T.cd8 = subset(pre.T,scibet.celltype %in% CD8subset)
saveRDS(pre.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/pre.T.cd8.rds')

post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/post.T.rds")
post.T.scibet = run.scibet(obj = post.T,ref = SciBet.ref)
saveRDS(post.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/post.T.scibet.rds')
post.T = AddMetaData(post.T,post.T.scibet)
post.T.cd8 = subset(post.T,scibet.celltype %in% CD8subset)
saveRDS(post.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_CTLA4_Moshe/post.T.cd8.rds')

#pd1---------------
pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/pre.T.rds")
pre.T.scibet = run.scibet(obj = pre.T,ref = SciBet.ref)
saveRDS(pre.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/pre.T.scibet.rds')
pre.T = AddMetaData(pre.T,pre.T.scibet)
pre.T.cd8 = subset(pre.T,scibet.celltype %in% CD8subset)
saveRDS(pre.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/pre.T.cd8.rds')

post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/post.T.rds")
post.T.scibet = run.scibet(obj = post.T,ref = SciBet.ref)
saveRDS(post.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/post.T.scibet.rds')
post.T = AddMetaData(post.T,post.T.scibet)
post.T.cd8 = subset(post.T,scibet.celltype %in% CD8subset)
saveRDS(post.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Moshe/post.T.cd8.rds')

# Melanoma_livnat---------------
post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/post.T.rds")
post.T.scibet = run.scibet(obj = post.T,ref = SciBet.ref)
saveRDS(post.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/post.T.scibet.rds')
post.T = AddMetaData(post.T,post.T.scibet)
post.T.cd8 = subset(post.T,scibet.celltype %in% CD8subset)
saveRDS(post.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/post.T.cd8.rds')

naive.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/naive.T.rds")
naive.T.scibet = run.scibet(obj = naive.T,ref = SciBet.ref)
saveRDS(naive.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/naive.T.scibet.rds')
naive.T = AddMetaData(naive.T,naive.T.scibet)
naive.T.cd8 = subset(naive.T,scibet.celltype %in% CD8subset)
saveRDS(naive.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/Melanoma_PD1_Livnat/naive.T.cd8.rds')

#nscle_liu--------------

pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/pre.T.rds")
pre.T.scibet = run.scibet(obj = pre.T,ref = SciBet.ref)
saveRDS(pre.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/pre.T.scibet.rds')
pre.T = AddMetaData(pre.T,pre.T.scibet)
pre.T.cd8 = subset(pre.T,scibet.celltype %in% CD8subset)
saveRDS(pre.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/pre.T.cd8.rds')

post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/post.T.rds")
post.T.scibet = run.scibet(obj = post.T,ref = SciBet.ref)
saveRDS(post.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/post.T.scibet.rds')
post.T = AddMetaData(post.T,post.T.scibet)
post.T.cd8 = subset(post.T,scibet.celltype %in% CD8subset)
saveRDS(post.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/post.T.cd8.rds')

#nscle_caushi_til
list <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Caushi/list.T.rds")
t.scibet.list = list()
seu.t.cd8.list = list()
for (i in 1:length(list)) {
  seu.t = list[[i]]
  sample.n = seu.t$sample[1]
  T.scibet = run.scibet(obj = seu.t,ref = SciBet.ref)
  t.scibet.list[[sample.n]] = T.scibet
  seu.t = AddMetaData(seu.t,T.scibet)
  seu.t.cd8 = subset(seu.t,scibet.celltype %in% CD8subset)
  seu.t.cd8.list[[sample.n]] = seu.t.cd8
}
saveRDS(t.scibet.list,file = "~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Caushi/scibet.list.rds")
saveRDS(seu.t.cd8.list,file = "~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Caushi/T.cd8.list.rds")

# scc---------

pre.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/pre.T.rds")
pre.T.scibet = run.scibet(obj = pre.T,ref = SciBet.ref)
saveRDS(pre.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/pre.T.scibet.rds')
pre.T = AddMetaData(pre.T,pre.T.scibet)
pre.T.cd8 = subset(pre.T,scibet.celltype %in% CD8subset)
saveRDS(pre.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/pre.T.cd8.rds')

post.T <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/post.T.rds")
post.T.scibet = run.scibet(obj = post.T,ref = SciBet.ref)
saveRDS(post.T.scibet,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/post.T.scibet.rds')
post.T = AddMetaData(post.T,post.T.scibet)
post.T.cd8 = subset(post.T,scibet.celltype %in% CD8subset)
saveRDS(post.T.cd8,file = '~/PanCancerICI/CD8TMetabolism/Data/CD8T/SCC_PD1/post.T.cd8.rds')

