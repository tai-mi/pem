# library -----------------------------------------------------------------

library(magrittr)
library(tidyverse)


# import ------------------------------------------------------------------

cl <- readRDS('data_misc/clinical.rds')


# load --------------------------------------------------------

require(immunarch)

rep1 <- repLoad('data_tils/')
rep1$meta %<>% mutate(patient=str_extract(Sample,'\\d+') %>% as.numeric()) %>%
  left_join(cl,by=c('patient'='sample')) %>% select(-Sample)
rep1$data %<>% set_names(rep1$meta$patient %>% as.character())

# filter all coding and inframe (doesn't change anything for these)
rep1$data <- map(rep1$data,~coding(.x) %>% inframes())

saveRDS(rep1,'data_misc/rep1.rds')


# exploration -------------------------------------------------------------

rep1 <- readRDS('data_misc/rep1.rds')
# total number of unique clonotypes per sample
repExplore(rep1$data,'volume','aa') %>% vis()+scale_y_log10()
# frequency distribution of clonotypes
repExplore(rep1$data,'count','aa') %>% vis()
# lengths of cdr3 seqs, not terribly useful
repExplore(rep1$data,'len','aa') %>% vis()+scale_y_log10()

repClonality(rep1$data,'clonal.prop',.perc=10) %>% vis()
repClonality(rep1$data,'top',.head=c(1,10,100)) %>% vis()

repOverlap(rep1$data,'overlap') %>% vis('heatmap2')
repOverlap(rep1$data,'jaccard') %>% vis('heatmap2')
# repOverlap(rep1$data,'tversky',.a=0.1,.b=0.5) %>% vis('heatmap2')
morRep <- repOverlap(rep1$data,'morisita')
morRep %>% vis('heatmap2')
# kinda a mess, jaccard/tversky are dominated by 14-18, overlap by 19-14, and morisita by 13-7

repOverlapAnalysis(morRep,'tsne') %>% vis()
repOverlapAnalysis(repOverlap(rep1$data,'jaccard'),'tsne') %>% vis()
repOverlapAnalysis(repOverlap(rep1$data,'overlap'),'tsne') %>% vis()

repDiversity(rep1$data,'inv.simp') %>% mutate(Sample=as.factor(Sample)) %>% ggplot(aes(Sample,Value,fill=Sample))+geom_col()
repSample(rep1$data,'downsample') %>% repDiversity('raref') %>% vis()
map(rep1$data,function(x){
  s2 <- 1000
  n <- ifelse(nrow(x)>s2,s2,nrow(x))
  x[sample(1:nrow(x),n,replace=F),]
}) %>% 
  repDiversity('raref',.step=1) %>% vis()
repDiversity(rep1$data,'hill') %>% vis()
repDiversity(rep1$data[-15],'hill') %>% vis()
repDiversity(rep1$data[-15],'hill',.min.q=4) %>% vis()

trackClonotypes(rep1$data,list(1,10)) %>% vis()
