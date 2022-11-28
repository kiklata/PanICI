library(reticulate)

conda_install("Renv", "scipy")
conda_install("Renv", "anndata==0.6.19")
conda_install("Renv", "loompy==2.0.17")
#conda_install("Renv", "scanpy")

use_condaenv("Renv")

# import SciPy (will use "r-reticulate" as per call to use_condaenv)

library(Seurat)
library(sceasy)

input_file <- 'anndata.h5ad.scanpy.hvg2000_PC10_res2.h5ad'
res_file = 'anndata.h5ad.scanpy.hvg2000_PC10_res2.rds'

scipy <- import("scipy")
anndata <- reticulate::import("anndata")
loompy <- reticulate::import('loompy')
#scanpy <- reticulate::import('scanpy')

scanpy = import('scanpy')

## Step 1: read input 
combined_seu = anndata$read_h5ad(input_file)


## Step 2: format transform and output 


# Seurat to AnnData
# sceasy::convertFormat(combined_seu, from="seurat", to="anndata", outFile=res_file)

# AnnData to Seurat
sceasy::convertFormat(obj = input_file, from="anndata", to="seurat", outFile =res_file)



