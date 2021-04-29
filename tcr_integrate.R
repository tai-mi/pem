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


# integrate to immunarch --------------------------------------------------

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

# scatterplot of proportion of each clone in one timepoint vs other
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


# expansion annotations ---------------------------------------------------

# Setup: run libraries, root, load rep2
rep2 <- readRDS('data_misc/rep2.rds')

# filter out patients without til data
withTils <- rep2$meta %>% filter(timepoint=='til') %>% .$patient
rep2$data <- rep2$data[rep2$meta$patient %in% withTils]
rep2$meta %<>% filter(patient %in% withTils)

map_dbl(rep2$data,nrow) %>% sort() 
#20-5 only has 2 doublets, 7-til has 1 8, 8-til has 2 8s, 7 8s, 3 8s
# might remove 20-5 from this analysis
# rep2$data <- rep2$data[names(rep2$data)!='20-5']
# rep2$meta %<>% filter(!((patient==20)&(timepoint==5)))

# annotate
dat@cell_tbl$expTil <- 0

for(p in withTils){
  for(i in c(1,2,8)){
    seqs <- rep2$data[[paste0(p,'-til')]] %>% 
      filter(Clones>=i) %>% .$CDR3.aa
    for(tp in rep2$meta %>% filter(patient==p,timepoint!='til') %>% .$timepoint){
      bars <- dat@contig_tbl %>% 
        filter(cdr3 %in% seqs,timepoint==tp,patient==p) %>% 
        .$barcode %>% unique()
      dat@cell_tbl$expTil[((dat@cell_tbl$barcode %in% bars)&
                             (dat@cell_tbl$patient==p)&
                             (dat@cell_tbl$timepoint==tp))] <- i
    }}}
# unnecessarilty assigns 1 level to 2's and 8's too
# could probably map2 somewhere
# saveRDS(dat,'data_misc/dat.rds')

# quant per til
map_dfr(withTils,function(p){
  temp1 <- map_dbl(c(1,2,8),function(i){
    rep2$data[[paste0(p,'-til')]] %>% 
      filter(Clones>=i) %>% nrow()
  })
  temp1[1] <- temp1[1]-temp1[2]
  temp1[2] <- temp1[2]-temp1[3]
  data.frame('patient'=p,'singlet'=temp1[1],'mid'=temp1[2],'expanded'=temp1[3])
}) %>% mutate(total=singlet+mid+expanded) %>% 
  mutate(across(c(expanded,mid,singlet),~.x/total)) %>% 
  pivot_longer(c('singlet','mid','expanded')) %>%
  arrange(patient) %>% mutate(patient=as.factor(patient)) %>%
  filter(name=='expanded') %>% #comment out to show all
  ggplot(aes(patient,value,fill=name))+geom_col(position='dodge') #position stack or fill to have stacked bars

