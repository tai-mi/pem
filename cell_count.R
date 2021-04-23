# library -----------------------------------------------------------------

library(magrittr)
library(tidyverse)


# load -------------------------------------------------------------------

stuff <- "35105,32291,31213,27492,25060,22267,20052,18650,18456,17122,16507,12946,10796,10134,9133,8758,8704,8591,4967,4725,3524,3336,3247,2966,2965,1706,1639,1248,1228,814,523,319,275,173,162,118" %>% 
  str_split(',') %>% unlist() %>% as.numeric()

typeFine <- c('mono14','t','nk','t','t','t','t','t','b','t','t','mono16','t',
              't','t','monoint','t','t','j','cdc1','mono','nk','nk','platelets',
              't','platelets','pdc','b','plasma','j','mono','mono','dc','j','j',
              'erythroid')
typeCoarse <- c('mono','t','nk','t','t','t','t','t','b','t','t','mono','t',
                't','t','mono','t','t','j','dc','mono','nk','nk','j',
                't','j','dc','b','b','j','mono','mono','dc','j','j',
                'erythroid')


# fail1 ----------------------------------------------------------------

cellNum <- tibble(type=typeFine,type2=typeCoarse,count=stuff)
cellNum %>% group_by(type) %>% summarise(n=sum(count))
# cellNum %>% group_by(type2) %>% summarise(n=sum(count)) %>% 
#   ggplot()+geom_col(aes(y=n,x=0,fill=type2),position='stack')
cellNum %>% filter(!(type2 %in% c('j','erythroid'))) %>% 
  group_by(type2) %>% summarise(n=sum(count)) %>% 
  ggplot()+geom_col(aes(x=1,y=n,fill=type2))+
  coord_polar('y',start=0)


# sourced -----------------------------------------------------------------

counts <- readRDS('data_sc/cell_source.rds') %>% 
  transmute(pem=p,healthy=h+s) %>% 
  mutate(type2=typeCoarse)
counts %>% pivot_longer(cols=c(pem,healthy)) %>% 
  filter(!(type2 %in% c('j','erythroid'))) %>%
  group_by(type2,name) %>% summarise(n=sum(value)) %>% 
  group_by(name) %>% summarise(n=n/sum(n),type2=type2) %>% 
  ggplot()+geom_col(aes(x=1,y=n,fill=type2))+
  coord_polar('y',start=0)+
  facet_wrap(vars(name))

# load patient level data
counts <- readRDS('data_sc/cell_orig.rds')
nm <- names(counts) %>% str_extract('(?<=^PEM)\\d+') %>% as.numeric() %>% 
  tibble(patient=.)
cl <- readRDS('data_misc/clinical.rds')
nm <- left_join(nm,cl,by=c('patient'='sample'))

# patient 24 doesn't have clinical data?
drop24 <- which(!is.na(nm$response))
nm <- nm[drop24,]
counts <- counts[,drop24]

piePlot <- function(df){
  df <- tibble(value=rowSums(df),type2=typeCoarse) %>% 
    group_by(type2) %>% summarize(value=sum(value)) %>% 
    filter(!(type2 %in% c('j','erythroid')))
  ggplot(df)+geom_col(aes(x=1,y=value,fill=type2))+
    coord_polar('y',start=0)
}

piePlot(counts[nm$lynch])
piePlot(counts[!nm$lynch])
piePlot(counts[nm$response %in% c('CR','PR')])
piePlot(counts[nm$response %in% c('PD','SD')])
piePlot(counts[nm$response=='CR'])
piePlot(counts[nm$response=='PR'])
piePlot(counts[nm$response=='SD'])
piePlot(counts[nm$response=='PD'])
