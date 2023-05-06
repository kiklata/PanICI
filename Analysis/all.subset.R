# get Compass Analysis DATA
All <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/All.rds")

# remove ICI NonTreat or ICI.efficacy NoInfo pts
all = subset(All,sample.timepoint != 'NonTreat')
all = subset(all,treatment.type != 'NonTreat')
all = subset(all,treatment.efficacy !='NoInfo')

# remove TPM dataset--Melanoma_Feldman
all = subset(all,DateType == 'Count')

# random select tumor samples in DATASET: NSCLC_Caushi for each patient
# max cell number

caushi.sample = names(table(all$sample.ID))[63:108]

selected.sample = c('NSCLC_CaushiMD01-004_tumor_1','NSCLC_CaushiMD01-005_tumor_7','NSCLC_CaushiMD01-010_tumor_1',
                    'NSCLC_CaushiMD01-019_tumor_1','NSCLC_CaushiMD01-024_tumor_1','NSCLC_CaushiMD043-003_tumor_3',
                    'NSCLC_CaushiMD043-006_tumor_1','NSCLC_CaushiMD043-008_tumor_1','NSCLC_CaushiMD043-011_tumor_3',
                    'NSCLC_CaushiNY016-007_tumor_1','NSCLC_CaushiNY016-014_tumor_1','NSCLC_CaushiNY016-015_tumor_2',
                    'NSCLC_CaushiNY016-021_tumor_1','NSCLC_CaushiNY016-022_tumor_1','NSCLC_CaushiNY016-025_tumor_2')

remove.sample = setdiff(caushi.sample,selected.sample)
all.sample = names(table(all$sample.ID))
sample.left = setdiff(all.sample,remove.sample)

all = subset(all,sample.ID %in% sample.left)

dim(all)

grep("^MT-", rownames(all),value = T)
grep("^HLA-", rownames(all),value = T)
grep('^RP[SL]', rownames(all),value = T)
grep('^MRP[SL]', rownames(all),value = T)

saveRDS(all,"~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/All.info.rds")

# convert geneSymbol using HGNChelper

library(HGNChelper)
exp.gene = rownames(All.info)

hgnc.check = checkGeneSymbols(exp.gene)
trans.gene = filter(hgnc.check,Approved == 'FALSE')
trans.gene = filter(trans.gene,Suggested.Symbol != '')

previous.exist.gene = trans.gene$x[trans.gene$Suggested.Symbol %in% exp.gene]
previous.non.gene = trans.gene$x[!trans.gene$Suggested.Symbol %in% exp.gene]

previous.exist.gene.new = trans.gene$Suggested.Symbol[trans.gene$Suggested.Symbol %in% exp.gene]
previous.non.gene.new = trans.gene$Suggested.Symbol[!trans.gene$Suggested.Symbol %in% exp.gene]

# previous.exist.gene.new have duplicated symbol, multi old symbol to one symbol
all.count = All.info@assays$RNA@counts

for (i in  1:length(previous.exist.gene)) {
  
  value.1 = all.count[previous.exist.gene[i],]
  value.2 = all.count[previous.exist.gene.new[i],]
  value = value.1 + value.2
  names(value) = NULL
  
  all.count = all.count[rownames(all.count) != previous.exist.gene[i],]
  all.count = all.count[rownames(all.count) != previous.exist.gene.new[i],]
  
  rrr = matrix(value,nrow = 1,ncol = length(value),dimnames = list(previous.exist.gene.new[i],names(value)))
  
  all.count = rbind(all.count,rrr)
  
  print(i)
}
saveRDS(all.count,file = 'all.count.remove1.rds')

# previous.non.gene check multi to one gene
all.count1 = all.count

for(i in 1:length(previous.non.gene)){
  
  value.1 = all.count1[previous.non.gene[i],]
  if (previous.non.gene.new[i] %in% rownames(all.count1)){
    value.2 = all.count1[previous.non.gene.new[i],]
  }else{
    value.2 = 0
  }
  value = value.1 + value.2
  
  all.count1 = all.count1[rownames(all.count1) != previous.non.gene[i],]
  
  if (previous.non.gene.new[i] %in% rownames(all.count1)){
  all.count1 = all.count1[rownames(all.count1) != previous.non.gene.new[i],]
  }
  
  rrr = matrix(value,nrow = 1,ncol = length(value),dimnames = list(previous.non.gene.new[i],names(value)))
  all.count1 = rbind(all.count1,rrr)
  print(i)
  
}
saveRDS(all.count1,file = 'all.count.symbol.correct.rds')

ALL.info <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/ALL.info.rds")

# hgnc

library(HGNChelper)
exp.gene = rownames(ALL.info)

hgnc.check = checkGeneSymbols(exp.gene)
left.gene = filter(hgnc.check,Approved == 'TRUE')$x

ALL.info = ALL.info[left.gene,]

saveRDS(ALL.info,file = 'ALL_CD8/ALL.info.hgnc.rds')

NSCLC_Liu.CD8 <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/NSCLC_PD1_Liu/NSCLC_Liu.CD8.rds")
ALL.info.hgnc <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/ALL.info.hgnc.rds")

metabolic.gene.name <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/metabolic.gene.name.rds")

NSCLC_Liu.CD8 = subset(NSCLC_Liu.CD8,treatment.type == 'PD1')

table(NSCLC_Liu.CD8$sample.ID)

# random select sample
sample.liu = names(table(NSCLC_Liu.CD8$sample.ID))
sample.liu = sample.liu[c(-1,-3,-8)]

NSCLC_Liu.CD8 = subset(NSCLC_Liu.CD8,sample.ID %in% sample.liu)

NSCLC_Liu.CD8$timepoint = NULL
NSCLC_Liu.CD8$CancerType = 'NSCLC'
NSCLC_Liu.CD8$Study = 'NSCLC_Liu'
NSCLC_Liu.CD8$Platform = '10X'
NSCLC_Liu.CD8$DateType = 'Count'
NSCLC_Liu.CD8$orig.ident = NULL

ALL.info.hgnc = merge(ALL.info.hgnc,NSCLC_Liu.CD8)
saveRDS(ALL.info.hgnc,file = 'ALL_CD8/ALL.info.hgnc.rds')

# metabolic gene

kegg = GSEABase::getGmt('~/PanCancerICI/CD8TMetabolism/Data/KEGG_metabolism_nc.gmt')

kegg.list = list()
for (i in 1:length(kegg)) {
  kegg.list[[kegg[[i]]@setName]] =  kegg[[i]]@geneIds
}

kegg.gene = kegg.list[[1]]

for (i in 2:length(kegg.list)) {
  gene = kegg.list[[i]]
  kegg.gene = append(kegg.gene,gene)
}

kegg.gene = kegg.gene[!duplicated(kegg.gene)]

#ALL.info.hgnc <- readRDS("~/PanCancerICI/CD8TMetabolism/Data/CD8T/ALL_CD8/ALL.info.hgnc.rds")
all.kegg = ALL.info.hgnc[kegg.gene,]
saveRDS(all.kegg,file = 'ALL_CD8/ALL.info.kegg.rds')

# shaozm
shaozmGene <- read_excel("~/PanCancerICI/CD8TMetabolism/Data/shaozmGene.xlsx")
colnames(shaozmGene) = shaozmGene[1,]
shaozmGene = shaozmGene[-1,]

# remove pseudogene by gtf------------------

gencode.v42.primary_assembly.annotation.gtf <- 
  read.delim2("~/bioinfo/reference/gencode.v42.primary_assembly.annotation.gtf.gz", 
              header=FALSE, comment.char="#")

gencode.v42.primary_assembly.annotation.gtf = 
  filter(gencode.v42.primary_assembly.annotation.gtf,V3 == 'gene')

gencode.v42.primary_assembly.annotation.gtf = gencode.v42.primary_assembly.annotation.gtf$V9
gencode.v42.primary_assembly.annotation.gtf = as.data.frame(gencode.v42.primary_assembly.annotation.gtf)
gtf = gencode.v42.primary_assembly.annotation.gtf

gtf.list = strsplit(gtf$gencode.v42.primary_assembly.annotation.gtf,";")

gtf.table.list = list()
for (i in 1:length(gtf.list)) {
  gtf.table = data.frame(geneid = gtf.list[[i]][1], 
                         gene_type = gtf.list[[i]][2], 
                         gene_name = gtf.list[[i]][3])
  gtf.table.list[[i]] = gtf.table
}
gtf.table = data.table::rbindlist(gtf.table.list)

gtf.table$geneid = substring(gtf.table$geneid,9)
gtf.table$gene_type = substring(gtf.table$gene_type,12)
gtf.table$gene_name = substring(gtf.table$gene_name,12)

# remove version number of gene ID
gtf.table$geneid = substring(gtf.table$geneid,1,15)

table(gtf.table$gene_type)

pseudogene = c('IG_pseudogene','IG_C_pseudogene','IG_J_pseudogene','IG_V_pseudogene','TR_V_pseudogene','TR_J_pseudogene',
               'Mt_tRNA_pseudogene','tRNA_pseudogene','snoRNA_pseudogene','snRNA_pseudogene',
               'scRNA_pseudogene','rRNA_pseudogene','misc_RNA_pseudogene','miRNA_pseudogene',
               'processed_pseudogene','pseudogene','transcribed_processed_pseudogene',
               'transcribed_unitary_pseudogene','transcribed_unprocessed_pseudogene',
               'translated_processed_pseudogene','translated_unprocessed_pseudogene',
               'unitary_pseudogene','unprocessed_pseudogene')
artifact = 'artifact'
TEC = 'TEC'

remove.gene.list = c(pseudogene,artifact,TEC)
all.gene.list = names(table(gtf.table$gene_type))
save.gene.list = setdiff(x = all.gene.list,y = remove.gene.list)
save.gene.name = filter(gtf.table, gene_type %in% save.gene.list)$gene_name

