
plot.sure.metastate = function(obj,cluster){
  library(dplyr)
  ptexpan = as.data.frame(table(obj$sample.ID,obj$MetaState))
  ptexpan = tidyr::spread(ptexpan,key = 'Var2',value = 'Freq')
  rownames(ptexpan) = ptexpan$Var1
  ptexpan$Var1 = NULL
  ptexpan$total = apply(ptexpan,1,sum)
  
  for (h in 1:(ncol(ptexpan)-1)) {
    ptexpan[,h] = ptexpan[,h]/ptexpan$total
  }
  
  ptefficacy = as.data.frame(table(obj$sample.ID,obj$treatment.efficacy))
  ptefficacy = tidyr::spread(ptefficacy,key = 'Var2',value = 'Freq')
  rownames(ptefficacy) = ptefficacy$Var1
  ptefficacy$Var1 = NULL
  ptefficacy$type = if_else(ptefficacy$NR==0,'R','NR')
  
  ptexpan$type = ptefficacy$type
  
  long.pt.pre <- tidyr::gather(ptexpan,key=MetaState,value = prop,-c('total','type'))
  long.pt.pre = filter(long.pt.pre,MetaState %in% cluster)
  
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  
  long.pt.pre$prop = long.pt.pre$prop + 0.00001
  long.pt.pre$type = factor(long.pt.pre$type, levels = c('R','NR'))
  p1 = ggboxplot(long.pt.pre, x = "MetaState", y = "prop",
                 fill = "type", palette = "nejm")+
    stat_compare_means(aes(group = type),
                       method = "wilcox.test",
                       label = "p.signif",
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "")))
  return(p1)
}
