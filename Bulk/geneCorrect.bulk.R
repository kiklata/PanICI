

gene.correct = function(obj) {
  library(dplyr)
  library(HGNChelper)
  
  exp.gene = rownames(obj)
  
  hgnc.check = checkGeneSymbols(exp.gene, species = 'human')
  trans.gene = filter(hgnc.check, Approved == 'FALSE')
  
  if (nrow(trans.gene) == 0) {
    return(obj)
  } else if (nrow(trans.gene) > 0) {
    trans.gene = filter(trans.gene, Suggested.Symbol != '')
    
    for (i in 1:nrow(trans.gene)) {
      trans.gene$Suggested.Symbol[i] = strsplit(trans.gene$Suggested.Symbol[i], ' /// ')[[1]][1]
      
    }
    
    previous.exist.gene = trans.gene$x[trans.gene$Suggested.Symbol %in% exp.gene]
    previous.non.gene = trans.gene$x[!trans.gene$Suggested.Symbol %in% exp.gene]
    
    if ((length(previous.exist.gene) == 0)) {
      previous.non.gene.new = trans.gene$Suggested.Symbol[!trans.gene$Suggested.Symbol %in% exp.gene]
      print('nonexist')
      
      pb2 = utils::txtProgressBar(style = 3)
      
      for (i in 1:length(previous.non.gene)) {
        value.1 = obj[previous.non.gene[i], ]
        if (previous.non.gene.new[i] %in% rownames(obj)) {
          value.2 = obj[previous.non.gene.new[i], ]
          value =  apply(value.1, 2, as.numeric) + apply(value.2, 2, as.numeric)
          
        } else{
          value.2 = 0
          value =  apply(value.1, 2, as.numeric)
          
        }
        
        obj = obj[rownames(obj) != previous.non.gene[i], ]
        
        if (previous.non.gene.new[i] %in% rownames(obj)) {
          obj = obj[rownames(obj) != previous.non.gene.new[i], ]
        }
        
        rrr = matrix(
          value,
          nrow = 1,
          ncol = length(value),
          dimnames = list(previous.non.gene.new[i], names(value))
        )
        colnames(rrr) = colnames(obj)
        
        obj = rbind(obj, rrr)
        
        utils::setTxtProgressBar(pb2, i / length(previous.non.gene))
      }
      close(pb2)
      
    } else if (length(previous.non.gene) == 0) {
      previous.exist.gene.new = trans.gene$Suggested.Symbol[trans.gene$Suggested.Symbol %in% exp.gene]
      print('exist')
      
      pb1 = utils::txtProgressBar(style = 3)
      
      for (i in  1:length(previous.exist.gene)) {
        value.1 = obj[previous.exist.gene[i], ]
        value.2 = obj[previous.exist.gene.new[i], ]
        value =  apply(value.1, 2, as.numeric) + apply(value.2, 2, as.numeric)
        
        obj = obj[rownames(obj) != previous.exist.gene[i], ]
        obj = obj[rownames(obj) != previous.exist.gene.new[i], ]
        
        rrr = matrix(
          value,
          nrow = 1,
          ncol = length(value),
          dimnames = list(previous.exist.gene.new[i], names(value))
        )
        obj = rbind(obj, rrr)
        
        utils::setTxtProgressBar(pb1, i / length(previous.exist.gene))
      }
      close(pb1)
      
    } else {
      previous.exist.gene.new = trans.gene$Suggested.Symbol[trans.gene$Suggested.Symbol %in% exp.gene]
      previous.non.gene.new = trans.gene$Suggested.Symbol[!trans.gene$Suggested.Symbol %in% exp.gene]
      
      
      # previous.exist.gene.new have duplicated symbol, multi old symbol to one symbol
      print('exist')
      
      pb1 = utils::txtProgressBar(style = 3)
      
      for (i in  1:length(previous.exist.gene)) {
        value.1 = obj[previous.exist.gene[i], ]
        value.2 = obj[previous.exist.gene.new[i], ]
        value =  apply(value.1, 2, as.numeric) + apply(value.2, 2, as.numeric)
        
        obj = obj[rownames(obj) != previous.exist.gene[i], ]
        obj = obj[rownames(obj) != previous.exist.gene.new[i], ]
        
        rrr = matrix(
          value,
          nrow = 1,
          ncol = length(value),
          dimnames = list(previous.exist.gene.new[i], names(value))
        )
        obj = rbind(obj, rrr)
        
        utils::setTxtProgressBar(pb1, i / length(previous.exist.gene))
      }
      close(pb1)
      # previous.non.gene check multi to one gene
      print('nonexist')
      
      pb2 = utils::txtProgressBar(style = 3)
      
      for (i in 1:length(previous.non.gene)) {
        value.1 = obj[previous.non.gene[i], ]
        if (previous.non.gene.new[i] %in% rownames(obj)) {
          value.2 = obj[previous.non.gene.new[i], ]
          value =  apply(value.1, 2, as.numeric) + apply(value.2, 2, as.numeric)
          
        } else{
          value.2 = 0
          value =  apply(value.1, 2, as.numeric)
          
        }
        
        obj = obj[rownames(obj) != previous.non.gene[i], ]
        
        if (previous.non.gene.new[i] %in% rownames(obj)) {
          obj = obj[rownames(obj) != previous.non.gene.new[i], ]
        }
        
        rrr = matrix(
          value,
          nrow = 1,
          ncol = length(value),
          dimnames = list(previous.non.gene.new[i], names(value))
        )
        colnames(rrr) = colnames(obj)
        
        obj = rbind(obj, rrr)
        
        utils::setTxtProgressBar(pb2, i / length(previous.non.gene))
      }
      close(pb2)
    }
    exp.gene = rownames(obj)
    
    hgnc.check = checkGeneSymbols(exp.gene)
    left.gene = filter(hgnc.check, Approved == 'TRUE')$x
    
    obj = obj[left.gene, ]
    
    return(obj)
  }
}
