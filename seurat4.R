# load main -----------------------------------------------------------------

library(dplyr)
library(purrr)
library(Matrix)
library(Seurat)
library(future)

options(future.globals.maxSize= 5*1024^4, future.rng.onMisuse="ignore") 

# registering URL/Wherever they are
setwd('/gpfs/ycga/scratch60/saltzman/hs685/santin/HHT')
file_list = list.dirs(full.names=FALSE, recursive = FALSE)
folder_list = paste0(list.dirs(full.names=TRUE,recursive=FALSE), "/outs/filtered_feature_bc_matrix/")

tmpList <- map( c(1:length(file_list)), function(tmp2){
  tmp <- Read10X(data.dir = folder_list[tmp2])
  tmp <- CreateSeuratObject(count=tmp, project=file_list[[tmp2]])
  tmp[['file_name']] <- file_list[[tmp2]]
  return(tmp)
})

tmpList <- merge(tmpList[[1]], y=tmpList[2:length(tmpList)] )
saveRDS(tmpList,'~/scratch60/pem/main.rds')


# load eric -------------------------------------------------------------

library(dplyr)
library(purrr)
library(Matrix)
library(future)
library(patchwork)
library(Seurat)

options(future.globals.maxSize= 5*1024^4, future.rng.onMisuse="ignore") 

# registering URL/Wherever they are
setwd('~/scratch60/pem/data/healthy')
folder_list = c('seurat10','seurat11','seurat12')
names_list = c('eric1','eric2','eric3')

tmpList <- map( c(1:length(folder_list)), function(tmp2){
  tmp <- Read10X(data.dir = folder_list[tmp2])
  tmp <- CreateSeuratObject(count=tmp, project=names_list[tmp2])
  tmp[['file_name']] <- folder_list[tmp2]
  return(tmp)
})

tmpList <- merge(tmpList[[1]], y=tmpList[2:length(tmpList)] )
saveRDS(tmpList,'~/scratch60/pem/eric.rds')


# load hd -----------------------------------------------------------------

library(dplyr)
library(purrr)
library(Matrix)
library(Seurat)
library(future)

options(future.globals.maxSize= 5*1024^4, future.rng.onMisuse="ignore") 

# registering URL/Wherever they are
setwd('~/scratch60/pem/data/healthy')
folder_list = c('HA5876BLD','HA5877BLD','HA5894BLD','HA5952BLD','HA5953BLD',
                'HA5957BLD')
names_list = c('HD1','HD2','HD3','HD4','HD5','HD6')

tmpList <- map( c(1:length(folder_list)), function(tmp2){
  tmp <- Read10X(data.dir = folder_list[tmp2])
  tmp <- CreateSeuratObject(count=tmp, project=names_list[tmp2])
  tmp[['file_name']] <- folder_list[tmp2]
  return(tmp)
})

tmpList <- merge(tmpList[[1]], y=tmpList[2:length(tmpList)] )
saveRDS(tmpList,'~/scratch60/pem/hd.rds')

# preprocess --------------------------------------------------------------

library(dplyr)
library(purrr)
library(Matrix)
library(patchwork)
library(Seurat)
library(future)
options(future.globals.maxSize= 5*1024^4, future.rng.onMisuse="ignore")
setwd('~/scratch60/pem/')

ldat <- list(
  readRDS('main.rds'),
  readRDS('eric.rds'),
  readRDS('hd.rds')
)

# ggplot(dat@meta.data)+geom_density(aes(nCount_RNA))+scale_x_log10()+geom_vline(xintercept = 400)
# ggplot(dat@meta.data)+geom_density(aes(nFeature_RNA))+scale_x_log10()+geom_vline(xintercept = 150)
# ggplot(hdat@meta.data)+geom_density(aes(propMt))+geom_vline(xintercept = 0.2)
# after checking, these look like good cutoffs, main has way more high mt cells (~25%) and has a bunch of cells around 500 nCount whereas the other two seem to have the cells with <500 nCount removed; both eric and main have nFeat much closer to 150 while hd is much higher peak

ldat <- map(ldat,function(dat){
  dat$propMt <- PercentageFeatureSet(dat,'^MT-')/100
  dat <- subset(dat, subset=(nCount_RNA>400 & nFeature_RNA>150 & propMt<0.2))
  dat <- SCTransform(dat,vars.to.regress=c('nCount_RNA','nFeatures_RNA','propMt'))
  dat
})
saveRDS(ldat,'ldat1.rds')

anchorFeats <- SelectIntegrationFeatures(ldat,nfeatures=3000)
ldat <- PrepSCTIntegration(ldat,'SCT',anchor.features=anchorFeats)
anchorNames <- FindIntegrationAnchors(ldat,'SCT',anchor.features=anchorFeats)
rm(anchorFeats)
rm(anchorNames)
ldat <- IntegrateData(anchorNames,normalization.method='SCT')
saveRDS(ldat,'ldat2.rds')

# integration -------------------------------------------------------------

library(dplyr)
library(purrr)
library(Matrix)
library(Seurat)
library(future)
options(future.globals.maxSize= 5*1024^4, future.rng.onMisuse="ignore")
setwd('~/scratch60/pem/')

# integrateSamples <- function(objList, nfeats=3000) {
#   anchorFeats <- SelectIntegrationFeatures(objList, nfeatures=nfeats)
#   objList %<>% PrepSCTIntegration(anchor.features=anchorFeats)
#   anchorNames <- FindIntegrationAnchors(objList,anchor.features=anchorFeats,normalization.method='SCT')
#   IntegrateData(anchorNames,normalization.method='SCT')
# }

ldat <- list(
  readRDS('main.rds'),
  readRDS('eric.rds'),
  readRDS('hd.rds')
)

anchorFeats <- SelectIntegrationFeatures(ldat,nfeatures=3000)
ldat <- PrepSCTIntegration(ldat,'SCT',anchor.features=anchorFeats)
anchorNames <- FindIntegrationAnchors(ldat,'SCT',anchor.features=anchorFeats)
rm(anchorFeats)
rm(anchorNames)
ldat <- IntegrateData(anchorNames,normalization.method='SCT')


# get data ----------------------------------------------------------------

library(dplyr)
library(purrr)
library(Matrix)
library(Seurat)
library(future)
options(future.globals.maxSize= 5*1024^4, future.rng.onMisuse="ignore")
setwd('~/scratch60/pem/')

ldat <- readRDS('ldat2.rds')
ldat <- RunPCA(ldat)
ldat <- FindNeighbors(ldat,dims=1:20)
saveRDS(ldat,'ldat3.rds')
ldat <- FindClusters(ldat,resolution=0.9)
ldat <- RunUMAP(ldat,dims=1:20)
saveRDS(ldat,'ldat4.rds')