# library -----------------------------------------------------------------

library(magrittr)
library(CellaRepertorium)
library(immunarch)
library(tidyverse)

rep3 <- readRDS('data_misc/rep3.rds')
cella <- readRDS('data_misc/cella3.rds')
cl <- readRDS('data_misc/clinical.rds')
dat <- readRDS('data_sc/meta_tdat.rds') %>% 
  mutate(patient=case_when(
    str_starts(orig.ident,'P')~str_extract(orig.ident,'(?<=PEM)\\d+'),
    str_starts(orig.ident,'eric')~orig.ident,
    orig.ident=='HD1'~'HA5876',orig.ident=='HD2'~'HA5877',
    orig.ident=='HD3'~'HA5894',orig.ident=='HD4'~'HA5952',
    orig.ident=='HD5'~'HA5953',orig.ident=='HD6'~'HA5957',
    T~NA_character_))


# exploratory -------------------------------------------------------------

# til overlap timepoints
repOverlap(rep3$data[rep3$meta$timepoint %in% c('0','7')],'jaccard') %>% 
  vis('heatmap2')
  # the pre-post tils always strong cluster in overlap/jaccard/morisita
  # 26/9/3/12 consistent, 19 and 7 iffy (7 really low counts)

# pbmc overlap timepoints
repOverlap(rep3$data[rep3$meta$timepoint %in% c('1','3','5')],'morisita') %>% 
  vis('heatmap2')
  # morisita looks beautiful, all patients cluster really strongly
  # some mixing with the low 20's in jaccard and overlap, but overall very strong

# overall overlap
repOverlap(rep3$data,'overlap') %>% vis('heatmap2')
  # overlap is nearly perfect (just 3til and 3pbmc separate but still have an island of overlap out there), jaccard good but mixes some up, morisita segregates like all the til from nontil

# change in overlap with pre/post tils over time
h1 <- repOverlap(rep3$data,'jaccard') %>% as.data.frame() %>% rownames_to_column('temp') %>% pivot_longer(-temp) %>% separate(temp,c('p1','t1'),'-') %>% separate(name,c('p2','t2'),'-') %>% na.omit() %>% filter(p1==p2,t1 %in% c('0','7'),t2 %in% c('1','3','5'))
ggplot(filter(h1,t1=='0'),aes(t2,value))+geom_boxplot(outlier.alpha=0)+geom_jitter(width=0.05,height=0)
ggplot(filter(h1,t1=='7'),aes(t2,value))+geom_boxplot(outlier.alpha=0)+geom_jitter(width=0.05,height=0)

# quantitative 1-5 change in til overlap
rbind(
  filter(h1,p1 %in% unique(h1$p1)[map_lgl(unique(h1$p1),~all(c('1','5') %in% filter(h1,p1==.x)$t2))],t1=='0') %>% pivot_wider(names_from=t2,values_from=value) %>% mutate(value=`5`-`1`),
  filter(h1,p1 %in% unique(h1$p1)[map_lgl(unique(h1$p1),~all(c('1','5') %in% filter(h1,p1==.x)$t2))],t1=='7') %>% pivot_wider(names_from=t2,values_from=value) %>% mutate(value=`5`-`1`)) %>% 
  ggplot(aes(t1,value,fill=t1))+geom_boxplot(outlier.alpha=0)+
  geom_jitter(width=0.05,height=0)
  # it kinda looked like something was there, but looks inconclusive

# cd8 cd4 stuff
h1 <- group_by(dat,patient,timepoint,clinical) %>% 
  summarize(cd8=mean(str_starts(cluster_coarse,'CD8')),
            cd4=mean(str_starts(cluster_coarse,'CD8',negate=T)))
ggplot(h1,aes(clinical,cd8,fill=clinical))+
  geom_boxplot(outlier.alpha=0)+
  geom_jitter(aes(shape=timepoint),width=0.05,height=0)
ggplot(h1,aes(clinical,cd8/cd4,fill=clinical))+
  geom_boxplot(outlier.alpha=0)+
  geom_jitter(aes(shape=timepoint),width=0.05,height=0)+
  scale_y_log10()
ggplot(filter(h1,timepoint!=''),aes(timepoint,cd8/cd4,fill=timepoint))+
  geom_boxplot(outlier.alpha=0)+
  geom_jitter(width=0.05,height=0)+
  scale_y_log10()
# no real overall trends

# cell types of clonotypes
  #' need to do cella, and then group patient-time-cdr3
  #' ID most freq cell type
  #' ID second most freq but just do this for top clones
nfilter <- 3
h1 <- add_clinical(cella@contig_tbl) %>% filter(clinical=='Healthy') %>% 
  group_split(patient,timepoint,chain,cdr3)
h1 <- h1[map_dbl(h1,~nrow(.x))>=nfilter]
h1 <- map_chr(h1,function(x){
  y <- na.omit(x$cluster_coarse) %>% as.character()
  if(length(y)<nfilter) NA_character_
  else{
    table1 <- table(y)
    if(length(table1)==1){
      sample(names(which(table1==max(table1))),1) %>% 
        paste0('-')
    } else{
      primary <- sample(names(which(table1==max(table1))),1)
      table1 <- table(y,exclude=primary)
      secondary <- sample(names(which(table1==max(table1))),1)
      paste0(primary,'-',secondary)
    }
  }
}) %>% na.omit() %>% 
  data.frame('temp'=.) %>% 
  separate(temp,c('primary','secondary'),'-') %>% 
  filter(primary!='Dead') %>% 
  mutate(secondary=case_when(secondary==''~'None',secondary=='Dead'~'None',T~secondary))
table(h1) %>% t() %>% as.data.frame.matrix() %>% 
  mutate(across(.fns=proportions)) %>% 
  rownames_to_column('secondary') %>% pivot_longer(-secondary) %>%
  rename(primary=name) %>% 
  # this is to replace the ones with cd8 naive as primary
  # mutate(secondary2=if_else(primary=='CD8 Naive','CD8 Naive',secondary),
  #        primary2=if_else(primary=='CD8 Naive',secondary,primary)) %>% 
  filter(primary!='CD8 Naive') %>% 
  ggplot(aes(primary,value,fill=secondary))+geom_col(position='stack')
  # don't have scaling worked out yet, and this is overall
  # but treg lot of overlap cd4act cd8, cd8exh-act strong, cd4 act-exh?
  # lot more treg and cd4 act in nlr, low treg in l, 

# lazy number of types per clonotype by clinical
add_clinical <- function(df,patient_col='patient'){
  p1 <- df[[patient_col]]
  df$clinical <- case_when(
    p1=='healthy'~'Healthy',
    p1 %in% c('14','25','5','2','6','23')~'Lynch-like responder',
    p1 %in% c('18','19','22','20','12','11','13')~'Nonlynch-like responder',
    p1 %in% c('1','3','7','8','9','10','15','16','17','21','24','26')~'Nonresponder',
    stringr::str_starts(p1,'(HD)|(HA)|(eric)')~'Healthy',
    T~NA_character_)
  if(any(is.na(df$clinical))) warning('invalid patient code')
  df
}
cella@contig_tbl %>% filter(chain=='TRB') %>% 
  group_by(patient,timepoint,cdr3) %>%
  summarize(ntypes=if_else(length(cluster_coarse)==1,NA_integer_,
                           length(unique(cluster_coarse)))) %>% 
  na.omit() %>% add_clinical() %>% 
  mutate(ntypes=factor(ntypes)) %>% 
  group_by(patient,timepoint,clinical,ntypes) %>% 
  summarize(count=length(ntypes)) %>% ungroup() %>% 
  group_by(patient,timepoint,clinical) %>% 
  summarize(prop=proportions(count),ntypes=ntypes) %>% ungroup() %>% 
  ggplot(aes(ntypes,prop))+geom_boxplot(aes(fill=clinical))+scale_y_sqrt()
  # not much clinical diff, maybe a bit more single in nr