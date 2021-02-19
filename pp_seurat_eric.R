# HEADER
library(batchtools)
library(data.table)
library(drake)
library(dplyr)
library(ggplot2)
library(hash)
library(listenv)
library(Matrix)
library(patchwork)
library(plotly)
library(Seurat)
library(tidyverse)

### FUTURE PARALLELIZATIONS X CPUs; future for Seurat framework and dopar for Ruddle archi
library(batchtools)
library(foreach)
library(future)
library(future.apply)
library(doFuture)
reg = makeRegistry(NA)
makeClusterFunctionsInteractive()
#library(future.batchtools
options(future.globals.maxSize= 5000000000000, future.rng.onMisuse="ignore") 

# registering URL/Wherever they are
file_list = list.dirs(path = "~/scratch60/jean/sub", full.names=FALSE, recursive = FALSE)
folder_list = paste0(list.dirs(path = "~/scratch60/jean/sub", full.names=TRUE,recursive=FALSE), "/filtered_feature_bc_matrix")

list_seurat_objects <- lapply( c(1:length(file_list)), function(tmp2){
    tmp <- Read10X(data.dir = folder_list[tmp2])
    tmp <- CreateSeuratObject(count=tmp, project=file_list[[tmp2]])
    tmp[['file_name']]<- file_list[[tmp2]]
    return(tmp)
})

####
big_seurat_object <- merge(list_seurat_objects[[1]], y=list_seurat_objects[2:length(list_seurat_objects)] )

big_seurat_object <- NormalizeData(big_seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(big_seurat_object)
big_seurat_object <- ScaleData(big_seurat_object, features = all.genes)

big_seurat_object <- FindVariableFeatures(object = big_seurat_object)
big_seurat_object <- RunPCA(big_seurat_object, features = VariableFeatures(object = big_seurat_object))

big_seurat_object <- FindNeighbors(big_seurat_object, dims = 1:20)
big_seurat_object <- FindClusters(big_seurat_object, resolution = 0.5)

big_seurat_object <- RunUMAP(big_seurat_object, dims = 1:20)

saveRDS(big_seurat_object, file = "~/scratch60/jean/big_seurat_object.rds")

