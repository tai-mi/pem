# libraries ---------------------------------------------------------------

library(magrittr)
library(tidyverse)


# main --------------------------------------------------------------------

h1 <- read_tsv('data_misc/combined_study_segments.seg')
h1 %<>% filter(chrom=='19')
h2 <- 55146317 #middle of LILRB1
filter(h1,loc.start<h2,loc.end>h2) %>% 
  ggplot(aes(reorder(ID,-seg.mean),seg.mean,fill=seg.mean>0))+geom_col()
