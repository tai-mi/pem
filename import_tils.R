# library -----------------------------------------------------------------

library(magrittr)
library(tidyverse)


# import ------------------------------------------------------------------

cl <- readRDS('data_misc/clinical.rds')


# load --------------------------------------------------------

require(immunarch)

rep1 <- repLoad('data_tils/') #misparses vj gene ties, will ignore for now
rep1$meta %<>% mutate(orig.file=Sample,
                      Sample=str_remove_all(Sample,'_'),
                      patient=as.numeric(str_extract(Sample,'(?<=PEM)\\d+')) %>% 
                        as.character(),
                      timepoint=if_else(str_detect(Sample,'R(?=TIL)'),'7','0'),
                      #timepoint sets pre as 0 and post as 7
                      group='til') %>%
  select(-Sample) 
rep1$data %<>% set_names(paste0(rep1$meta$patient,'-',rep1$meta$timepoint))

# filter all coding and inframe 
rep1$data <- map(rep1$data,~coding(.x) %>% inframes())
# map_dbl(rep1$data,~nrow(.x)) %>% sort() #both for 7 are low, kinda 8-0 and 3-0
# map_dbl(rep1$data,~sum(.x$Clones)) %>% sort()

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
h1 <- repOverlap(rep1$data,'morisita')
h1 %>% vis('heatmap2')
# strong cluster of 0-7 pairs for all of them; with 7s removed,a few pairs, but not clinically related

repOverlapAnalysis(h1,'tsne') %>% vis()
repOverlapAnalysis(repOverlap(rep1$data,'jaccard'),'tsne') %>% vis()
repOverlapAnalysis(repOverlap(rep1$data,'overlap'),'tsne') %>% vis()

repDiversity(rep1$data,'inv.simp') %>% mutate(Sample=as.factor(Sample)) %>% 
  ggplot(aes(Sample,Value,fill=Sample))+geom_col()+scale_y_log10()+
  theme(axis.text.x=element_text(angle=90),legend.pos='none')
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
