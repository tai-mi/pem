# library -----------------------------------------------------------------

library(magrittr)
library(tidyverse)


# import ------------------------------------------------------------------

cl <- readRDS('clinical.rds')


# immunarch method --------------------------------------------------------

require(immunarch)

rep1 <- repLoad('data_tils/')
rep1$meta %<>% mutate(patient=str_extract(Sample,'\\d+') %>% as.numeric()) %>% 
  left_join(cl,by=c('patient'='sample')) %>% select(-Sample)
rep1$data %<>% set_names(rep1$meta$patient %>% as.character())

# filter all coding and inframe (doesn't change anything for these)
rep1$data <- map(rep1$data,~coding(.x) %>% inframes())

# total number of unique clonotypes per sample
repExplore(rep1$data,'volume','aa') %>% vis()+scale_y_log10()
# frequency distribution of clonotypes
repExplore(rep1$data,'count','aa') %>% vis()
# lengths of cdr3 seqs, not terribly useful
repExplore(rep1$data,'len','aa') %>% vis()+scale_y_log10()


