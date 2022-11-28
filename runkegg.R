# run.kegg----------------------

run.kegg = function(obj){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  genetoanno = obj$gene
  genelist.ENTREZID = bitr(genetoanno,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,2]
  
  ekegg = enrichKEGG(genelist.ENTREZID,organism = 'hsa',keyType = "kegg",
                     pvalueCutoff = 0.05,pAdjustMethod = "BH")
  return(ekegg)
}
