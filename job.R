require(tidyverse)

spade <- dat@contig_tbl %>% group_split(patient,timepoint) %>% 
  map(~.x$cdr3 %>% table() %>% unlist() %>% 
        SpadeR::Diversity('abundance'))