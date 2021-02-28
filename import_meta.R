
# libraries ---------------------------------------------------------------

library(tidyverse)


# import ------------------------------------------------------------------

meta <- readRDS('data_sc/metadata.rds') %>% 
  as_tibble() %>% 
  select(origIdent=orig.ident,ncount=nCount_RNA,nfeat=nFeature_RNA,
         cluster=seurat_clusters,s=S.Score,g2m=G2M.Score,phase=Phase,propMt)

ggplot(meta)+geom_density(aes(x=ncount))+scale_x_log10()
ggplot(meta)+geom_density(aes(x=nfeat))+scale_x_log10()
ggplot(meta)+geom_density(aes(x=g2m))
ggplot(meta)+geom_density(aes(x=s))
ggplot(meta)+geom_density(aes(x=propMt))
