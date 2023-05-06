heatmap(scale(table(Before.tsne.kegg$scibet.celltype,Before.tsne.kegg$RNA_snn_res.0.4)),Rowv = NA,Colv = NA)
DimPlot(Before.tsne.kegg,group.by = 'RNA_snn_res.0.4')
all.time.plot[["Before"]][["0.4"]]