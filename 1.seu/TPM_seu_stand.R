all.T.raw <- readRDS("~/PaperCD8/data/unannotated/all.T.raw.rds")

tpm = subset(all.T.raw,DateType == 'TPM')

tpm[["percent.mt"]] <- PercentageFeatureSet(tpm, pattern = "^MT-")
tpm = subset(tpm,percent.mt<20) # 

source("~/PaperCD8/code/geneblacklist.R")
rpgene = grep('^RP[LS]',rownames(tpm),value = T)
mtgene = 'MALAT1'
final.gene = setdiff(rownames(tpm),y = c(immgene,rpgene,mtgene))
tpm = tpm[final.gene,]

# remove HSP gene
HSPgene = grep('^HSP',rownames(tpm),value = T)
final.gene = setdiff(rownames(tpm),HSPgene)
tpm = tpm[final.gene,]

# filter <3 gene, <200 cell
count = tpm[['RNA']]@counts
seu = CreateSeuratObject(counts = count,min.cells = 3,min.features = 100,meta.data = tpm@meta.data)
saveRDS(seu,file = 'tpm.rds')

tpm <- readRDS("~/PaperCD8/data/tpm.rds")
data = log(tpm@assays$RNA@counts+1)
data = Matrix::as.matrix(data,Class = 'dgcMatrix')
tpm <- SetAssayData(object = tpm, slot = "data", new.data = data)

tpm = tpm %>% FindVariableFeatures() %>% ScaleData()

tpm = tpm %>% RunPCA() %>% RunHarmony(group.by.vars = 'sample.ID',assay.use='RNA', 
                                      plot_convergence = F,theta = 2.5, 
                  kmeans_init_nstart=20, kmeans_init_iter_max=5000) %>% 
RunUMAP(reduction = 'harmony',dims = 1:30,reduction.name = 'umap_harmony') %>% 
FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
FindClusters(resolution = seq(0.2,1,0.2))

saveRDS(tpm,file = 'tpm.harmony.rds')
