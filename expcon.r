# libraries ---------------------------------------------------------------

source('figures_source.r')
load_tcr()
dat <- load_main('tc')


# uhhhhhhhhhh -------------------------------------------------------------

#' so like um yeah what's the goal here
#' are we going framework or specific questions
#' and each of those is gonna necessitate some of the other
#' what are our questions that we specifically want to look at
#' TM expansion in NLR, change in peripheral clonality lynch
#' yeah so framework
#' we need to know what our measurement even is
#' categorical is easiest but also a lie
#' definitely should be downsampling but that's prob bad for a lot of them
#' ok so no downsampling or maybe on an individual patient level but that's also not ideal since diff patients having diff numbers that they're downsampling to could be confounding and slight towards one or the other but also they're already slanted w/o downsampling
#' what if instead of downsample just scale counts by the diff in total counts between the two giving one of the pair fractional counts
#' nah i'm just not gonna downsample - all of this is single cell so it should be fine?
#' soooo we doing categories or what - singlet, 2-7, 8+ (this should really be scaled by the total number of cells but i'm ignoring that)
#' could also do 1-2 vs 3+
#' using this framework, gets a network of transitions - focus on the most expanded from one and the other and how they transfer

# full join of clonotypes with classification of pre and post expansion category
expcon_compare_categorical <- function(df1,df2,overlap8=F){
  h1 <- full_join(
    mutate(df1,clones1=case_when(Clones==1~'1',Clones %in% 2:7~'2',Clones>7~'8')) %>% 
      select(cdr3=CDR3.aa,clones1),
    mutate(df2,clones2=case_when(Clones==1~'1',Clones %in% 2:7~'2',Clones>7~'8')) %>% 
      select(cdr3=CDR3.aa,clones2),
    by='cdr3'
  ) %>% mutate(across(everything(),~tidyr::replace_na(.,'0')),
               expcon=interaction(clones1,clones2))
  if(overlap8){
    int8 <- sum(h1$expcon=='8.8')
    h1 <- c(int8/sum(h1$clones1=='8'),int8/sum(h1$clones2=='8'))
  } 
  h1
}
# expcon_compare_categorical(rep4$data$`14-1`,rep4$data$`14-5`) %>% 
#   select(clones1,clones2) %>% table()

# calc proportional change in df2 of the top n clonotypes from df1
expcon_compare_rank <- function(df1,df2,n=10){
  left_join(
    arrange(df1,desc(Clones)) %>% mutate(proportion1=proportions(Clones)) %>% 
      slice(1:n) %>% select(cdr3=CDR3.aa,proportion1),
    mutate(df2,proportion2=proportions(Clones)) %>% 
      select(cdr3=CDR3.aa,proportion2),
    by='cdr3'
  ) %>% mutate(across(everything(),~tidyr::replace_na(.,0)),
               expcon=(proportion2-proportion1)/proportion1)
}

with15 <- group_by(rep4$meta,patient) %>% 
  summarize(a=if_else(all(c('1','5') %in% timepoint),patient[1],NA_character_)) %>%
  na.omit() %>% .$patient
# rep4$data[(rep4$meta$patient %in% with15)&(rep4$meta$timepoint %in% c(1,5))] %>% map_dbl(nrow) %>% sort() # 23 and 20 are really low, 5-1 is 459 and all others >500
with15 <- with15[with15 %!in% c('23','20')]
h1 <- map_dfr(with15,function(i){
  size1 <- 25
  h1 <- expcon_compare_rank(rep4$data[[paste0(i,'-1')]],
                            rep4$data[[paste0(i,'-5')]],n=size1)$expcon
  data.frame(patient=i,rank=1:size1,expcon=h1)
}) %>% add_clinical() %>% mutate(rank=as.factor(rank))
annot_row <- data.frame(patient=with15) %>% add_clinical() %>% 
  column_to_rownames('patient')
pivot_wider(h1,id_cols=c(patient,clinical),names_from=rank,values_from=expcon) %>% 
  arrange(clinical) %>% select(-clinical) %>% 
  column_to_rownames('patient') %>% 
  pheatmap::pheatmap(cluster_rows=F,cluster_cols=F,annotation_row=annot_row,
                     color=c(colorRampPalette(c('blue','white'))(100*abs(range(h1$expcon)[1])),colorRampPalette(c('white','red'))(100*abs(range(h1$expcon)[2]))))
# looks uhhh pretty fucking uninterpretable where the hell do i go with this?

# 1-5 overlap of expanded downsampled
map(with15,function(i){
  x <- list(rep4$data[[paste0(i,'-1')]],rep4$data[[paste0(i,'-5')]]) %>% 
    repSample()
  expcon_compare_categorical(x[[1]],x[[2]],overlap8=T)
}) %>% as.data.frame() %>% t() %>% 
  set_colnames(c('c1','c2')) %>% as_tibble() %>% 
  mutate(patient=with15) %>% add_clinical() %>% 
  plotBoxplot(value_col='c2')
# bootstrapped
h1 <- map_dfr(1:50,function(ii){
  map(with15,function(i){
    x <- list(rep4$data[[paste0(i,'-1')]],rep4$data[[paste0(i,'-5')]]) %>% 
      repSample()
    expcon_compare_categorical(x[[1]],x[[2]],overlap8=T)
  }) %>% as.data.frame() %>% t() %>% 
    set_colnames(c('c1','c2')) %>% as_tibble() %>% 
    mutate(patient=with15) %>% add_clinical() %>% 
    group_by(clinical) %>% summarize(v=mean(c2))
})
plotBoxplot(h1,'v')

# tm stuff
map(with15,function(i){
  x <- list(rep4$data[[paste0(i,'-1')]],rep4$data[[paste0(i,'-5')]]) %>% 
    map(~filter(.x,!is.na(expTil)))
  expcon_compare_categorical(x[[1]],x[[2]],overlap8=T)
}) %>% as.data.frame() %>% t() %>% 
  set_colnames(c('c1','c2')) %>% as_tibble() %>% 
  mutate(patient=with15) %>% add_clinical() %>% 
  plotBoxplot(value_col='c2')
# bootstrapped
h1 <- map_dfr(1:50,function(ii){
  map(with15,function(i){
    x <- list(rep4$data[[paste0(i,'-1')]],rep4$data[[paste0(i,'-5')]]) %>% 
      map(~filter(.x,!is.na(expTil))) %>% repSample()
    expcon_compare_categorical(x[[1]],x[[2]],overlap8=T)
  }) %>% as.data.frame() %>% t() %>% 
    set_colnames(c('c1','c2')) %>% as_tibble() %>% 
    mutate(patient=with15) %>% add_clinical() %>% 
    group_by(clinical) %>% summarize(v=mean(c2,na.rm=T))
})
plotBoxplot(h1,'v')


## new era
expcon_overlap_rank <- function(df1,df2,n=25){
  if((n>nrow(df1))|(n>nrow(df2))) warning('Clones less than n')
  cdr3.1 <- arrange(df1,desc(Clones))$CDR3.aa[1:min(n,nrow(df1))]
  cdr3.2 <- arrange(df2,desc(Clones))$CDR3.aa[1:min(n,nrow(df2))]
  sum(cdr3.1 %in% cdr3.2)/length(unique(c(cdr3.1,cdr3.2)))
}

# 1-5 overlap of expanded downsampled
map_dbl(with15,function(i){
  x <- list(rep4$data[[paste0(i,'-1')]],rep4$data[[paste0(i,'-5')]]) %>% 
    repSample()
  expcon_overlap_rank(x[[1]],x[[2]],n=25)
}) %>% data.frame('value'=.,patient=with15) %>% add_clinical() %>% 
  plotBoxplot(value_col='value')
# bootstrapped
h1 <- map_dfr(1:50,function(ii){
  map_dbl(with15,function(i){
    x <- list(rep4$data[[paste0(i,'-1')]],rep4$data[[paste0(i,'-5')]]) %>% 
      repSample()
    expcon_overlap_rank(x[[1]],x[[2]],n=25)
  }) %>% data.frame('value'=.,patient=with15) %>% add_clinical() %>% 
    group_by(clinical) %>% summarize(value=mean(value))
})
plotBoxplot(h1,'value')
# tm bootstrapped
h1 <- map_dfr(1:50,function(ii){
  map_dbl(with15,function(i){
    x <- list(rep4$data[[paste0(i,'-1')]],rep4$data[[paste0(i,'-5')]]) %>% 
      map(~filter(.x,!is.na(expTil))) %>% repSample()
    expcon_overlap_rank(x[[1]],x[[2]],n=15)
  }) %>% data.frame('value'=.,patient=with15) %>% add_clinical() %>% 
    group_by(clinical) %>% summarize(value=mean(value,na.rm=T))
})
plotBoxplot(h1,'value')
