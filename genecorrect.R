gene.correct = function(obj){
  
  library(dplyr)
  library(Seurat)
  library(HGNChelper)
  
  exp.gene = rownames(obj)
  
  hgnc.check = checkGeneSymbols(exp.gene,species = 'human')
  trans.gene = filter(hgnc.check,Approved == 'FALSE')
  trans.gene = filter(trans.gene,Suggested.Symbol != '')
  
  for (i in 1:nrow(trans.gene)) {
    
    trans.gene$Suggested.Symbol[i] = strsplit(trans.gene$Suggested.Symbol[i],' /// ')[[1]][1]
    
  }
  
  previous.exist.gene = trans.gene$x[trans.gene$Suggested.Symbol %in% exp.gene]
  previous.non.gene = trans.gene$x[!trans.gene$Suggested.Symbol %in% exp.gene]
  
  previous.exist.gene.new = trans.gene$Suggested.Symbol[trans.gene$Suggested.Symbol %in% exp.gene]
  previous.non.gene.new = trans.gene$Suggested.Symbol[!trans.gene$Suggested.Symbol %in% exp.gene]
  
  # previous.exist.gene.new have duplicated symbol, multi old symbol to one symbol
  all.count = obj@assays$RNA@counts
  print('exist')
  
  pb1 = utils::txtProgressBar(style = 3)

  for (i in  1:length(previous.exist.gene)) {
    
    
    
    value.1 = all.count[previous.exist.gene[i],]
    value.2 = all.count[previous.exist.gene.new[i],]
    value = value.1 + value.2
    names(value) = NULL
    
    all.count = all.count[rownames(all.count) != previous.exist.gene[i],]
    all.count = all.count[rownames(all.count) != previous.exist.gene.new[i],]
    
    rrr = matrix(value,nrow = 1,ncol = length(value),dimnames = list(previous.exist.gene.new[i],names(value)))
    
    all.count = rbind(all.count,rrr)
    
    utils::setTxtProgressBar(pb1,i/length(previous.exist.gene))
  }
  close(pb1)
  # previous.non.gene check multi to one gene
  all.count1 = all.count
  print('nonexist')
  
  pb2 = utils::txtProgressBar(style = 3)
  
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
    
    utils::setTxtProgressBar(pb2,i/length(previous.non.gene))
  }
  close(pb2)
  
  correct.obj = CreateSeuratObject(all.count1,meta.data = obj@meta.data)
  
  exp.gene = rownames(correct.obj)
  
  hgnc.check = checkGeneSymbols(exp.gene)
  left.gene = filter(hgnc.check,Approved == 'TRUE')$x

  correct.obj = correct.obj[left.gene,]
  
  return(correct.obj)

}
