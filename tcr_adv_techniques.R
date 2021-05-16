# libraries ---------------------------------------------------------------

library(magrittr)
library(immunarch)
library(CellaRepertorium)
library(tidyverse)


# load --------------------------------------------------------------------

rep <- readRDS('data_misc/rep2_healthy.rds')
dat <- readRDS('data_misc/dat_healthy.rds')

rep$meta %<>% mutate(coarse=str_replace(timepoint,'1|3|5','pbmc'))


# immunarch ---------------------------------------------------------------

# V gene overlap
trv <- geneUsage(rep$data,.gene='hs.trbv',.quant='count',.ambig='exc',
                  .norm=T)
# trv <- geneUsage(rep$data,.gene='hs.trav',.quant='count',.ambig='exc',
#                   .norm=T)

geneUsageAnalysis(trv,.method='tsne') %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  ggplot(aes(DimI,DimII,label=rowname,color=rep$meta$patient))+geom_point()+
  ggrepel::geom_text_repel()
geneUsageAnalysis(trv,.method='pca+hclust') %>% vis()
geneUsageAnalysis(trv,.method='js') %>% vis(.plot='heatmap2')
geneUsageAnalysis(trv,.method='js')[c(20:32,35:48,50:51,53:82),c(20:32,35:48,50:51,53:82)] %>% pheatmap::pheatmap()
geneUsageAnalysis(trv,.method='mds+dbscan')$data %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  ggplot(aes(DimI,DimII,label=rowname,color=rep$meta$coarse))+geom_point()+
  ggrepel::geom_text_repel()

# rep overlap
jc <- repOverlap(rep$data,'jaccard')
mor <- repOverlap(rep$data,'morisita')

vis(jc,.plot='heatmap2')
vis(mor,.plot='heatmap2')

repOverlapAnalysis(jc,'tsne') %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  ggplot(aes(DimI,DimII,label=rowname,color=rep$meta$coarse))+geom_point()+
  ggrepel::geom_text_repel()
repOverlapAnalysis(mor,'tsne') %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  ggplot(aes(DimI,DimII,label=rowname,color=rep$meta$lynch))+geom_point()+
  ggrepel::geom_text_repel()
repOverlapAnalysis(mor,'tsne+kmeans',.k=6) %>% vis()
#### try aggregating patients and then running this analysis on them as wholes

vdj <- data.table::fread('data_rep/vdjdb_single.tsv')
#### do stuff

# kmers
k5 <- getKmers(rep$data,5)
# k8 <- getKmers(rep$data,8)

vis(k5,.head=30,.position='fill')
kmer_profile(k5[[1]],.method='freq') %>% vis()
kmer_profile(k5[[1]],.method='freq') %>% vis(.plot='seq')
