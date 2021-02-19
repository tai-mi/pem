
# library -----------------------------------------------------------------

library(SingleR)
library(tidyverse)

# stuff -------------------------------------------------------------------

coarse <- readRDS('sc_data/singler_coarse.rds')
fine <- readRDS('sc_data/singler_fine.rds')

coarse$main <- coarse$labels %>% str_extract('.*(?= \\()')
coarse$fine <- coarse$labels %>% str_extract('(?<=\\().*(?=\\))')
coarse$labels <- c()
fine$main <- fine$labels %>% str_extract('.*(?= \\()')
fine$fine <- fine$labels %>% str_extract('(?<=\\().*(?=\\))')
fine$labels <- c()

countCoarse <- coarse$labels %>% table() %>% sort(decreasing=T)
countFine <- fine$labels %>% table() %>% sort(decreasing=T)

# shit, the coarse is random and the fine is all Bcell
