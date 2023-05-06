library(uwot)

umap.res = uwot::umap(score.cohort1.norm)

head(umap.res)

plot(umap.res,col=score.cohort1.norm$response,pch=16,asp = 1,
     xlab = "UMAP_1",ylab = "UMAP_2",
     main = "A UMAP visualization of the Metabolism ssGSEA scores")
abline(h=0,v=0,lty=2,col="gray")
legend("topright",title = "Response",inset = 0.01,
       legend = unique(score.cohort1.norm$response),pch=16,
       col = unique(score.cohort1.norm$response))


seu = CreateSeuratObject(t(score.cohort1.norm))
icb = subset(raw_cohort1_tumor,expansion !='n/a')
seu = seu[,colnames(icb)]
seu = AddMetaData(seu,icb@meta.data)

seu = subset(seu,timepoint =='Pre')

seu = ScaleData(seu)

seu = FindVariableFeatures(seu)
seu = FindNeighbors(seu)
seu = FindClusters(seu,resolution = 0.1)

seu = RunUMAP(seu,dims = 1:30)

DimPlot(seu)
DimPlot(seu,split.by = 'expansion')
table(seu$RNA_snn_res.0.1,seu$expansion)
allmarker = FindAllMarkers(seu)
