# library -----------------------------------------------------------------

library(magrittr)
library(CellaRepertorium)
library(immunarch)
library(tidyverse)

rep1 <- readRDS('data_misc/rep1.rds')
dat <- readRDS('data_misc/dat_healthy.rds')
cl <- readRDS('data_misc/clinical.rds')


# root --------------------------------------------------------------------

# takes df (already filtered to one set to analyze)
toImmunarch <- function(df){
  nTot <- nrow(df)
  df <- add_count(df,cdr3)
  df[!duplicated(df$cdr3),] %>% rename(
    Clones=n,
    CDR3.nt=cdr3_nt,
    CDR3.aa=cdr3,
    V.name=v_gene,
    D.name=d_gene,
    J.name=j_gene
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

combining <- function(dat1){
  dat1@cell_tbl$combined <- 
    paste(dat1@cell_tbl$barcode,
          dat1@cell_tbl$patient %>% str_replace('^[a-zA-Z].*','healthy'),
          tidyr::replace_na(dat1$cell_tbl$timepoint,''),
          sep='_')
  dat1@contig_tbl$combined <- 
    paste(dat1@contig_tbl$barcode,
          dat1@contig_tbl$patient %>% str_replace('^[a-zA-Z].*','healthy'),
          tidyr::replace_na(dat1$contig_tbl$timepoint,''),
          sep='_')
  dat1
}


# old integrate dat to immunarch -------------------------------------------

# make an immunarch of everything
dat2 <- dat@contig_tbl %>% group_split(patient,timepoint) %>% 
  map(~toImmunarch(.x) %>% mutate(timepoint=replace_na(timepoint,'h')))
names(dat2) <- map(dat2, function(x) c(x$patient[1],x$timepoint[1])) %>% 
  map(~paste(.x,collapse='-')) %>% unlist()
# set as patient-timepoin, healthy timepoint is h and patient name is distict

rep2 <- c(rep1$data,dat2)
rm(rep1)
meta2 <- names(rep2) %>% str_split('-') %>% as.data.frame() %>% t() %>% 
  as.data.frame() %>% set_colnames(c('patient','timepoint')) %>% 
  left_join(cl,by='patient')
rep2 <- list(meta2,rep2) %>% set_names(c('meta','data'))

rep2$data %<>% map(function(x){
  x$V.name[x$V.name=='None'] <- NA
  x$D.name[x$D.name=='None'] <- NA
  x$J.name[x$J.name=='None'] <- NA
  x
})

saveRDS(rep2,'data_misc/rep2.rds')


# expansion annotations ---------------------------------------------------

rep2 <- readRDS('data_misc/rep2.rds')

# filter out patients without til data
withTils <- rep2$meta %>% filter(timepoint=='0') %>% .$patient
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


# integrate healthy ---------------------------------------------------------


# Setup: run libraries, root, load rep2
rep2 <- readRDS('data_misc/rep2.rds')
cl <- readRDS('data_misc/clinical.rds') ### warning: cl has changed
cl$sample <- as.character(cl$sample)

# filter out patients without til data
withTils <- rep2$meta %>% filter(timepoint=='til') %>% .$patient

# convert dat_healthy to immunarch
dat2 <- dat@contig_tbl %>% group_split(patient,timepoint) %>% 
  map(~toImmunarch(.x))
names(dat2) <- map(dat2, function(x) c(x$patient[1],if_else(is.na(x$timepoint[1]),'h',x$timepoint[1]))) %>% 
  map(~paste(.x,collapse='-')) %>% unlist()

names(rep1$data) %<>% paste0('-til')
rep2 <- c(rep1$data,dat2)
meta2 <- names(rep2) %>% str_split('-') %>% as.data.frame() %>% t() %>% 
  as.data.frame() %>% set_colnames(c('patient','timepoint')) %>% 
  left_join(cl,by=c('patient'='sample'))
rep2 <- list(meta2,rep2) %>% set_names(c('meta','data'))

rep2$data %<>% map(function(x){
  x$V.name[x$V.name=='None'] <- NA
  x$D.name[x$D.name=='None'] <- NA
  x$J.name[x$J.name=='None'] <- NA
  x
})
saveRDS(rep2,'data_misc/rep2_healthy.rds')


# exp healthy full --------------------------------------------------------

dat <- readRDS('data_misc/dat_healthy.rds')
dat@cell_tbl$expTil <- 0
dat@cell_tbl$expPbmc <- 0

#### add data for healthy pbmc level
rep2 <- readRDS('data_misc/rep2_healthy.rds')
rep2$data <- rep2$data[rep2$meta$timepoint=='h']
rep2$meta %<>% filter(timepoint=='h')

map_dbl(rep2$data,nrow) %>% sort() 

for(p in rep2$meta$patient %>% unique()){
  for(i in c(1,2,8)){
    seqs <- rep2$data[[paste0(p,'-h')]] %>% 
      filter(Clones>=i) %>% .$CDR3.aa
    bars <- dat@contig_tbl %>% 
      filter(cdr3 %in% seqs,patient==p) %>% 
      .$barcode %>% unique()
    dat@cell_tbl$expPbmc[((dat@cell_tbl$barcode %in% bars)&
                           (dat@cell_tbl$patient==p))] <- i
    }}

#### add pem til data
rep2 <- readRDS('data_misc/rep2_healthy.rds')

withTils <- rep2$meta %>% filter(timepoint=='til') %>% .$patient
rep2$data <- rep2$data[rep2$meta$patient %in% withTils]
rep2$meta %<>% filter(patient %in% withTils)

map_dbl(rep2$data,nrow) %>% sort() 

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

#### add pem pbmc data
rep2 <- readRDS('data_misc/rep2_healthy.rds')

rep2$data <- rep2$data[rep2$meta$timepoint %in% c('1','3','5')]
rep2$meta %<>% filter(timepoint %in% c('1','3','5'))

map_dbl(rep2$data,nrow) %>% sort() 
# 23-1 and 20-5 have like nothing, 21-1 kinda low too (200)

for(p in rep2$meta$patient %>% unique()){
  for(i in c(1,2,8)){
    for(tp in rep2$meta %>% filter(patient==p) %>% .$timepoint){
      seqs <- rep2$data[[paste0(p,'-',tp)]] %>% 
        filter(Clones>=i) %>% .$CDR3.aa
      bars <- dat@contig_tbl %>% 
        filter(cdr3 %in% seqs,timepoint==tp,patient==p) %>% 
        .$barcode %>% unique()
      dat@cell_tbl$expPbmc[((dat@cell_tbl$barcode %in% bars)&
                             (dat@cell_tbl$patient==p)&
                             (dat@cell_tbl$timepoint==tp))] <- i
    }}}
saveRDS(dat,'data_misc/dat_healthy.rds')

# new integration 7.20 -----------------------------------------------------

#' uhhhhh so likeee what are we even doing here I mean I guess uhhhhhhhhhhhh
#' hm so let's focus we need to do
#' 190477/197854 cells (with contigs)
  #' many cells are still in cell_tbl but have contigs that have been filtered out
  #' 8528 duplicated (so at most 8.5k with two TRBs, could be some with 3+)

# setup
rep1 <- readRDS('data_misc/rep1.rds')
dat <- readRDS('data_misc/dat_healthy.rds')
dat@contig_tbl$expTil <- c()
dat@contig_tbl$expPbmc <- c()
dat@cell_tbl$expTil <- c()
dat@cell_tbl$expPbmc <- c()
dat@contig_tbl %<>% mutate(timepoint=replace_na(timepoint,'h'))
dat@cell_tbl %<>% mutate(timepoint=replace_na(timepoint,'h'))
# add expPbmc to dat (all clones get same number, ie: all 80 matching contigs have expPbmc=80)
dat@contig_tbl %<>% add_count(patient,timepoint,chain,cdr3,name='expPbmc')
# create df of seqs and expTil
h1 <- map_dfr(names(rep1$data),function(x){
  rep1$data[[x]] %>% mutate(patient=str_extract(x,'.*(?=-)'),
                            timepoint=str_extract(x,'.$'))
}) %>% add_count(patient,timepoint,CDR3.aa,wt=Clones,name='expTil') %>% 
  select(CDR3.aa,patient,timepoint,expTil) %>%
  unique() #add_count leaves the duplicated in, need to remove
# add expTil and expPost to pbmc
h1 <- left_join(filter(dat@contig_tbl,chain=='TRB'),
                filter(h1,timepoint=='0') %>% select(-timepoint),
                by=c('patient'='patient','cdr3'='CDR3.aa')) %>% 
  left_join(filter(h1,timepoint=='7') %>% rename(expPost=expTil) %>% 
              select(-timepoint),
            by=c('patient'='patient','cdr3'='CDR3.aa'))
dat@contig_tbl <- rbind(filter(dat@contig_tbl,chain=='TRA') %>% 
                          mutate(expTil=NA_integer_,expPost=NA_integer_),
                        h1)
rm(h1) #42018/199005 TRBs have expTil, 8018/42589 TRB have expPost
# map contig characteristics to cells
  # have 28634 duplicate contigs - up to that many cells with mult TRB
  # only mapping contig chars, sort by expPbmc more accurate to represent cdr3
dat <- canonicalize_cell(dat,contig_filter_args=chain=='TRB',
                         tie_break_keys=c('expPbmc','umis'),
                         contig_fields=c('cdr3','cdr3_nt','v_gene','d_gene',
                                         'j_gene','c_gene','expPbmc'))
  #190477 have cdr3 - all of them with TRB contigs
# add expansion annotations to cells
  # canonicalize doesn't necessarily retain highest exptil/post since sort by pbmc
dat <- canonicalize_cell(dat,contig_filter_args=chain=='TRB',
                         tie_break_keys='expTil',contig_fields='expTil')
dat <- canonicalize_cell(dat,contig_filter_args=chain=='TRB',
                         tie_break_keys='expPost',contig_fields='expPost')
#### first, branch this and commit, load rep1 into seurat from this branch
# add repdata to seurat
 #' copy the whole seurat matching scheme, and add all cell_tbl to seurat meta
 #' for this whole section, dat1 is seurat, dat is cella
dat1@meta.data <- left_join(
  mutate(dat1@meta.data,
         temp=str_remove(rownames(dat1@meta.data),'(?<=-1).*'),
         timepoint=tidyr::replace_na(timepoint,'h')),
  mutate(dat@cell_tbl,
         patient=case_when(str_starts(patient,'\\d')~patient,T~'healthy')) %>% 
    select(barcode,patient,timepoint,expPbmc,expTil,expPost,cdr3,cdr3_nt,
           v_gene,d_gene,j_gene,c_gene),
  by=c('temp'='barcode','patient'='patient','timepoint'='timepoint'))
dat@cell_tbl <- left_join(
  dat@cell_tbl,
  mutate(dat1@meta.data,
         temp=str_remove(colnames(dat1),'(?<=-1).*'),
         timepoint=tidyr::replace_na(timepoint,'h'),
         patient=case_when(
           str_starts(orig.ident,'PEM')~str_extract(orig.ident,'(?<=PEM)\\d+'),
           orig.ident=='HD1'~'HA5876',orig.ident=='HD2'~'HA5877',
           orig.ident=='HD3'~'HA5894',orig.ident=='HD4'~'HA5952',
           orig.ident=='HD5'~'HA5953',orig.ident=='HD6'~'HA5957',
           T~orig.ident)) %>% 
    select(temp,timepoint,patient,integrated_snn_res.0.5,
           integrated_snn_res.0.8,integrated_snn_res.1.1,
           cluster_fine,cluster_coarse),
  by=c('barcode'='temp','patient'='patient','timepoint'='timepoint'))
# transfer seurat from cell to contig
dat@contig_tbl %<>% left_join(
  select(dat@cell_tbl,barcode,patient,timepoint,integrated_snn_res.0.5,
         integrated_snn_res.0.8,integrated_snn_res.1.1,
         cluster_fine,cluster_coarse),
  by=c('barcode','patient','timepoint'))
# pseudo-bulk contigs and add to immunarch
choose1 <- function(x){#randomly chooses among most frequent of char_vec
  x <- na.omit(x)
  if(length(x)==0) NA
  else if(length(x)==1) x
  else{
    table1 <- table(x)
    sample(names(which(table1==max(table1))),1)
  }
}
rep2 <- dat@contig_tbl %>% filter(chain=='TRB') %>% #this is slow, sorry
  select(patient,timepoint,cdr3,v_gene,d_gene,j_gene,cdr3_nt,expPbmc,
         expTil,expPost,integrated_snn_res.0.5,integrated_snn_res.0.8,
         integrated_snn_res.1.1,cluster_coarse,cluster_fine) %>% 
  group_split(patient,timepoint,cdr3) %>% 
  map_dfr(function(y){
    if(nrow(y)==1){
      rename(y,V.name=v_gene,D.name=d_gene,J.name=j_gene,CDR3.nt=cdr3_nt)
    } else{
      data.frame(
        y[1,c('patient','timepoint','cdr3','expPbmc','expTil','expPost')],
        'V.name'=choose1(y$v_gene),'D.name'=choose1(y$d_gene),
        'J.name'=choose1(y$j_gene),'CDR3.nt'=choose1(y$cdr3_nt),
        'integrated_snn_res.0.5'=choose1(y$integrated_snn_res.0.5),
        'integrated_snn_res.0.8'=choose1(y$integrated_snn_res.0.8),
        'integrated_snn_res.1.1'=choose1(y$integrated_snn_res.1.1),
        'cluster_coarse'=choose1(y$cluster_coarse),
        'cluster_fine'=choose1(y$cluster_fine)
      )
    }
  }) %>% rename(Clones=expPbmc,CDR3.aa=cdr3) %>% 
  mutate(D.start=NA_integer_,D.end=NA_integer_,V.end=NA_integer_,
         J.start=NA_integer_,VJ.ins=NA_integer_,DJ.ins=NA_integer_,
         Sequence=NA_character_,
         V.name=if_else(V.name=='None',NA_character_,V.name),
         D.name=if_else(D.name=='None',NA_character_,D.name),
         J.name=if_else(J.name=='None',NA_character_,J.name)) %>% 
  group_split(patient,timepoint) %>%
  map(~mutate(.x,Proportion=proportions(Clones)))
names(rep2) <- map_chr(rep2,~paste0(.x$patient[1],'-',.x$timepoint[1]))
rep2 <- list('data'=c(rep1$data,rep2))
rep2$meta <- tibble(patient=str_extract(names(rep2$data),'.*(?=-)'),
                    timepoint=str_extract(names(rep2$data),'.$')) %>% 
  mutate(clinical=case_when(
    patient=='healthy'~'Healthy',
    patient %in% c('14','25','5','2','6','23')~'Lynch-like responder',
    patient %in% c(18,19,22,20,12,11,13)~'Nonlynch-like responder',
    T~'Nonresponder'
  ))
