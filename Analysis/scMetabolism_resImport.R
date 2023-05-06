normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
# BC
raw_cohort1_tumor <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/raw_cohort1_tumor.rds")
score.cohort1 <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/AntiPD1_BC/score.cohort1.rds")
raw_cohort1_tumor = subset(raw_cohort1_tumor,expansion !='n/a')
pre = subset(raw_cohort1_tumor,timepoint == 'Pre')
on = subset(raw_cohort1_tumor,timepoint =='On')
score.pre = score.cohort1[,colnames(pre)]
score.on = score.cohort1[,colnames(on)]

score.cohort1.pre.norm = apply(t(score.pre), 2, normalize)
score.cohort1.on.norm = apply(t(score.on), 2, normalize)

save(score.pre,score.on,score.cohort1.pre.norm,score.cohort1.on.norm,
     file = '~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/score.cohort1.rdata')


score.cohort2 <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/score.cohort2.rds")
raw_cohort2_tumor <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/raw_cohort2_tumor.rds")
pre = subset(raw_cohort2_tumor,timepoint == 'Pre')
on = subset(raw_cohort2_tumor,timepoint =='On')
score.pre = score.cohort2[,colnames(pre)]
score.on = score.cohort2[,colnames(on)]

score.cohort2.pre.norm = apply(t(score.pre), 2, normalize)
score.cohort2.on.norm = apply(t(score.on), 2, normalize)

save(score.pre,score.on,score.cohort2.pre.norm,score.cohort2.on.norm,
     file = '~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/score.cohort2.rdata')


score.cohort2 <- readRDS("AntiPD1_BC/score.cohort2.rds")
score.cohort2.norm = apply(t(score.cohort2), 2, normalize)
saveRDS(score.cohort2.norm,file = 'AntiPD1_BC/score.cohort2.norm.rds')

score.melanoma <- readRDS("GSE115978_Melanoma/score.melanoma.rds")
raw_tumor <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/GSE115978_Melanoma/raw_tumor.rds")
table(raw_tumor$treatment.group)
post = subset(raw_tumor,treatment.group =='post.treatment')
naive = subset(raw_tumor,treatment.group == 'treatment.naive')
score.post = score.melanoma[,colnames(post)]
score.naive = score.melanoma[,colnames(naive)]
score.melanoma.post.norm = apply(t(score.post), 2, normalize)
score.melanoma.naive.norm = apply(t(score.naive), 2, normalize)

save(score.melanoma.post.norm,score.melanoma.naive.norm,file = 'GSE115978_Melanoma/score.melanoma.rdata')


raw_BCC_tumor <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/GSE123813_BCC/raw_BCC_tumor.rds")
score.bcc <- readRDS("GSE123813_BCC/score.bcc.rds")
pre = subset(raw_BCC_tumor,treatment == 'pre')
post = subset(raw_BCC_tumor,treatment =='post')
score.pre = score.bcc[,colnames(pre)]
score.post = score.bcc[,colnames(post)]
score.pre.norm = apply(t(score.pre), 2, normalize)
score.post.norm = apply(t(score.post), 2, normalize)
save(score.pre,score.post,score.pre.norm,score.post.norm,
     file = '~/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/score.bcc.rdata')


raw_tumor <- readRDS("~/PanCancerICI/TumorMetabolism/Data/scData/SCP1288_ccRCC/raw_tumor.rds")
score.ccRCC <- readRDS("SCP1288_ccRCC/score.ccRCC.rds")

icb = subset(raw_tumor,ICB_Response %in% c('ICB_PD', 'ICB_PR', 'ICB_SD'))
score.ccRCC = score.ccRCC[,colnames(icb)]
score.ccRCC.norm = apply(t(score.ccRCC), 2, normalize)
saveRDS(score.ccRCC.norm,file = '~/PanCancerICI/TumorMetabolism/Data/scData/SCP1288_ccRCC/score.ccRCC.norm.rds')
