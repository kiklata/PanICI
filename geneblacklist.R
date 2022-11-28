#blacklist
library(dplyr)
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

#immunoglobin gene/TCR gene

immunoglobulin = c('IG_C_gene','IG_C_pseudogene','IG_D_gene','IG_J_gene','IG_J_pseudogene',
                   'IG_pseudogene','IG_V_gene','IG_V_pseudogene')

immgene = filter(gtf.table,gene_type %in% immunoglobulin)$gene_name



tcr = c('TR_C_gene','TR_D_gene','TR_J_gene',
        'TR_J_pseudogene', 'TR_V_gene','TR_V_pseudogene')

tcrgene = filter(gtf.table,gene_type %in% tcr)$gene_name

# RP gene
#rpgene = grep('^RP[LS]',rownames(ALL.info.hgnc),value = T)
