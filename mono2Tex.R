# monocle
library(monocle,lib.loc = '/home/shpc_100839/R/x86_64-pc-linux-gnu-library/4.2')
library(Seurat)
library(SeuratWrappers)
library(dplyr)

CD8.Tex.harmony <- readRDS("~/PaperCD8/data/Tex/CD8.downsample.rds")

expr_matrix = GetAssayData(CD8.Tex.harmony, assay = 'RNA', slot = 'counts')
pdata = CD8.Tex.harmony@meta.data
fdata = data.frame(gene_short_name=row.names(CD8.Tex.harmony), row.names = row.names(CD8.Tex.harmony))

pd = new('AnnotatedDataFrame', data = pdata)
fd = new('AnnotatedDataFrame', data = fdata)

cds <- newCellDataSet(as(expr_matrix, "sparseMatrix"),
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(cds[expressed_genes,],
                                      fullModelFormulaStr = "~manual.celltype.minor")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.05))
cds <- setOrderingFilter(cds, ordering_genes) 

## Trajectory step 2: reduce data dimensionality
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds) 

saveRDS(cds,file = 'CD8.cds.rds')
saveRDS(ordering_genes,file = 'CD8.genes.rds')
