# library -----------------------------------------------------------------

library(magrittr)
library(CellaRepertorium)
library(immunarch)
library(tidyverse)

rep1 <- readRDS('data_misc/rep1.rds')
dat <- readRDS('data_misc/dat.rds')


# root --------------------------------------------------------------------

# takes df (already filtered to one set to analyze)
toImmunarch <- function(df){
  nTot <- nrow(df)
  df <- add_count(df,cdr3)
  df[!duplicated(df$cdr3),] %>% select(
    Clones=n,
    CDR3.nt=cdr3_nt,
    CDR3.aa=cdr3,
    V.name=v_gene,
    D.name=d_gene,
    J.name=j_gene,
    patient,
    timepoint
  ) %>% mutate(
    V.end=NA,
    D.start=NA,
    D.end=NA,
    J.start=NA,
    VJ.ins=NA,
    VD.ins=NA,
    DJ.ins=NA,
    Sequence=NA,.after=6) %>% 
    mutate(Proportion=Clones/nTot,.after=1) %>% 
    arrange(desc(Clones))
}

# tracks clonotypes across pbmc and til timepoints
trackCl2 <- function(p,...){
  require(immunarch)
  p <- as.character(p)
  temp <- filter(dat@contig_tbl,patient==p) %>%
    group_split(timepoint) %>%
    map(~toImmunarch(.x))
  temp %<>% set_names(filter(dat@contig_tbl,patient==p) %>%
                        group_by(timepoint) %>% group_keys() %>% unlist())
  if(p %in% names(rep1$data)) temp$til <- rep1$data[[p]]
  trackClonotypes(temp,...)
}

# make an immunarch of everything
dat2 <- dat@contig_tbl %>% group_split(patient,timepoint) %>% 
  map(~toImmunarch(.x))
names(dat2) <- map(dat2, function(x) c(x$patient[1],x$timepoint[1])) %>% 
  map(~paste(.x,collapse='-')) %>% unlist()

names(rep1$data) %<>% paste0('-til')
rep2 <- c(rep1$data,dat2)
meta2 <- names(rep2) %>% str_split('-') %>% as.data.frame() %>% t() %>% 
  as.data.frame() %>% set_colnames(c('patient','timepoint')) %>% 
  mutate(patient=as.numeric(patient)) %>% 
  left_join(cl,by=c('patient'='sample'))
rep2 <- list(meta2,rep2) %>% set_names(c('meta','data'))

rep2$data %<>% map(function(x){
  x$V.name[x$V.name=='None'] <- NA
  x$D.name[x$D.name=='None'] <- NA
  x$J.name[x$J.name=='None'] <- NA
  x
})

saveRDS(rep2,'data_misc/rep2.rds')

# expcon -----------------------------------------------------------------

rep2 <- readRDS('data_misc/rep2.rds')

propComp <- function(p,x,y){
  # will error if patient is missing required data
  joined <- full_join(rep2$data[[paste0(p,'-',x)]],
                      rep2$data[[paste0(p,'-',y)]],
                      by=c('CDR3.aa'))
  # joined$Proportion.x[is.na(joined$Proportion.x)] <- 0
  joined$Proportion.y[is.na(joined$Proportion.y)] <- 0
  ggplot(joined,aes(Proportion.x,Proportion.y))+
    geom_jitter()+scale_x_log10()+scale_y_log10()+
    labs(title=paste(p,'-',x,'vs',y))
}

#' quantify proportion/prop of top clones from til that present in pbmc
#' correlations between them?
#' scatterplot of pbmc1 vs til