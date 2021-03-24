library(future)
library(Seurat)

options(future.globals.maxSize = 1500 * 1024^3)
setwd('~/scratch60/pem')

dat <- readRDS('dat_pca.rds')

dat <- RunUMAP(dat,dims=1:20)
saveRDS(dat,'dat_umap.rds')
UMAPPlot(dat)
ggplot2::ggsave('features/umap_sct.png',width=12,height=10,dpi=300,units='in')

dat <- FindNeighbors(dat,dims=1:20)
saveRDS(dat,'dat_neighbors.rds')

dat <- FindClusters(dat,resolution=0.7)
saveRDS(dat,'dat_clust.rds')