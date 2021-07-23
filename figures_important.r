# library1 -----------------------------------------------------------------

library(magrittr)
library(CellaRepertorium)
library(immunarch)
library(tidyverse)

rep3 <- readRDS('data_misc/rep3.rds')
cella <- readRDS('data_misc/cella3.rds')
cl <- readRDS('data_misc/clinical.rds')
dat <- readRDS('data_sc/meta_tdat.rds') %>% 
  mutate(patient=case_when(
    str_starts(orig.ident,'P')~str_extract(orig.ident,'(?<=PEM)\\d+'),
    T~orig.ident))


# sample overlap ----------------------------------------------------------

# repoverlap til only
repOverlap(rep3$data[rep3$meta$timepoint %in% c('0','7')],'jaccard') %>% 
  vis('heatmap2')
# repoverlap all tcr
pheatmap::pheatmap(
  repOverlap(rep3$data,'overlap'),
  color=colorRampPalette(c("#67001f","#d6604d","#f7f7f7",
                           "#4393c3","#053061"))(1000),
  breaks=c(seq(min(h1,na.rm=T),0.09, length.out=1001),max(h1,na.rm=T)),
  annotation_col=data.frame('clinical'=factor(rep3$meta$clinical),
                            row.names=names(rep3$data)))


# cell types of clonotypes ------------------------------------------------

#### modify to add partial counts when stuff tied for primary
nfilter <- 3 #only clones with at least this many IDed cells
h1 <- add_clinical(cella@contig_tbl) %>% filter(clinical=='Healthy') %>% 
  group_split(patient,timepoint,chain,cdr3)
h1 <- h1[map_dbl(h1,~nrow(.x))>=nfilter]
h1 <- map_chr(h1,function(x){
  y <- na.omit(x$cluster_coarse) %>% as.character()
  if(length(y)<nfilter) NA_character_
  else{
    table1 <- table(y)
    if(length(table1)==1){
      sample(names(which(table1==max(table1))),1) %>% 
        paste0('-')
    } else{
      primary <- sample(names(which(table1==max(table1))),1)
      table1 <- table(y,exclude=primary)
      secondary <- sample(names(which(table1==max(table1))),1)
      paste0(primary,'-',secondary)
    }
  }
}) %>% na.omit() %>% 
  data.frame('temp'=.) %>% 
  separate(temp,c('primary','secondary'),'-') %>% 
  filter(primary!='Dead') %>% 
  mutate(secondary=case_when(secondary==''~'None',secondary=='Dead'~'None',T~secondary))
table(h1) %>% t() %>% as.data.frame.matrix() %>% 
  mutate(across(.fns=proportions)) %>% 
  rownames_to_column('secondary') %>% pivot_longer(-secondary) %>%
  rename(primary=name) %>% 
  # this is to replace the ones with cd8 naive as primary
  # mutate(secondary2=if_else(primary=='CD8 Naive','CD8 Naive',secondary),
  #        primary2=if_else(primary=='CD8 Naive',secondary,primary)) %>% 
  # filter(primary!='CD8 Naive') %>% 
  ggplot(aes(primary,value,fill=secondary))+geom_col(position='stack')
# don't have scaling to correct for the number of cells of each type yet
# but treg lot of overlap cd4act cd8, cd8exh-act strong, cd4 act-exh?
# lot more treg and cd4 act in nlr, low treg in l, 
