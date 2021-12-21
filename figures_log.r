# library -----------------------------------------------------------------

library(magrittr)
library(CellaRepertorium)
library(immunarch)
library(tidyverse)

rep4 <- readRDS('data_misc/rep4.rds')
rep4$data <- map(rep4$data,~mutate(.x,Proportion=proportions(.x$Clones)))
cella <- readRDS('data_misc/cella4.rds')
cl <- readRDS('data_misc/clinical.rds')
# dat <- readRDS('data_sc/meta_ldat.rds') %>% 
dat <- readRDS('data_sc/meta_tc.rds') %>% 
  mutate(patient=case_when(
    str_starts(orig.ident,'P')~str_extract(orig.ident,'(?<=PEM)\\d+'),
    str_starts(orig.ident,'eric')~orig.ident,
    orig.ident=='HD1'~'HA5876',orig.ident=='HD2'~'HA5877',
    orig.ident=='HD3'~'HA5894',orig.ident=='HD4'~'HA5952',
    orig.ident=='HD5'~'HA5953',orig.ident=='HD6'~'HA5957',
    T~NA_character_))

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
cella@contig_tbl %<>% add_clinical()
cella@cell_tbl %<>% add_clinical()

plotBoxplot <- function(df,value_col='Value',patient_col1='patient'){
  # if(!('clinical' %in% names(df))) dat <- add_clinical(dat,patient_col1)
  df$Value <- df[[value_col]]
  ggplot(df,aes(clinical,Value,fill=clinical))+
    geom_boxplot(outlier.alpha=0)+geom_jitter(width=0.05,height=0)
}
require(ggprism)
cols_clinical <- c('Lynch-like responder'='#D25B3F','Nonlynch-like responder'='#356870','Nonresponder'='#8D7233','Healthy'='#3C367E')
cols_clinical_dark <- colorspace::darken(cols_clinical,0.6) %>% set_names(c('Lynch-like responder','Nonlynch-like responder','Nonresponder','Healthy'))
prism <- function(gg){
  gg+theme_prism()+scale_fill_manual(values=cols_clinical)+
    scale_color_manual(values=cols_clinical_dark)
}


# temporal overlap -----------------------------------------------------------

# til overlap timepoints
repOverlap(rep4$data[rep4$meta$timepoint %in% c('0','7')],'jaccard') %>% 
  vis('heatmap2')
  # the pre-post tils always strong cluster in overlap/jaccard/morisita
  # 26/9/3/12 consistent, 19 and 7 iffy (7 really low counts)

# pbmc overlap timepoints
repOverlap(rep4$data[rep4$meta$timepoint %in% c('1','3','5')],'morisita') %>% 
  vis('heatmap2')
  # morisita looks beautiful, all patients cluster really strongly
  # some mixing with the low 20's in jaccard and overlap, but overall very strong

# overall overlap
repOverlap(rep4$data,'overlap') %>% vis('heatmap2')
  # overlap is nearly perfect (just 3til and 3pbmc separate but still have an island of overlap out there), jaccard good but mixes some up, morisita segregates like all the til from nontil

# change in overlap with pre/post tils over time
h1 <- repOverlap(rep4$data,'jaccard') %>% as.data.frame() %>% rownames_to_column('temp') %>% pivot_longer(-temp) %>% separate(temp,c('p1','t1'),'-') %>% separate(name,c('p2','t2'),'-') %>% na.omit() %>% filter(p1==p2,t1 %in% c('0','7'),t2 %in% c('1','3','5'))
# replace 0 with 7 in filter to run post-treatment
ggplot(filter(h1,t1=='0'),aes(t2,value))+geom_boxplot(outlier.alpha=0)+geom_jitter(width=0.05,height=0)
filter(h1,t1=='0',t2!='3') %>% pivot_wider(id_cols=p1,names_from=t2,values_from=value) %>% 
  na.omit() %>% mutate(value=`5`-`1`) %>% # this is only 5-1 diff
  add_clinical('p1') %>% 
  ggplot(aes(reorder(p1,desc(value)),value,fill=clinical))+geom_col()
  # looks like generally stronger 5-7 (but small); tendency of responders to large changes 1-5 in 0 overlap but changes could be pos or neg (there's some paper on this?)

# quantitative 1-5 change in til overlap
rbind(
  filter(h1,p1 %in% unique(h1$p1)[map_lgl(unique(h1$p1),~all(c('1','5') %in% filter(h1,p1==.x)$t2))],t1=='0') %>% pivot_wider(names_from=t2,values_from=value) %>% mutate(value=`5`-`1`),
  filter(h1,p1 %in% unique(h1$p1)[map_lgl(unique(h1$p1),~all(c('1','5') %in% filter(h1,p1==.x)$t2))],t1=='7') %>% pivot_wider(names_from=t2,values_from=value) %>% mutate(value=`5`-`1`)) %>% 
  ggplot(aes(t1,value,fill=t1))+geom_boxplot(outlier.alpha=0)+
  geom_jitter(width=0.05,height=0)
  # it kinda looked like something was there, but looks inconclusive


# cd4-cd8 stuff -----------------------------------------------------------

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
pivot_wider(filter(h1,timepoint!=''),id_cols=patient,names_from=timepoint,values_from=cd8) %>% 
  transmute(patient=patient,change=`5`-`1`) %>% na.omit() %>% add_clinical() %>% 
  ggplot(aes(reorder(patient,desc(change)),change,fill=clinical))+geom_col()
# no real overall trends


# cell types of clonotypes ------------------------------------------------

nfilter <- 3
h1 <- add_clinical(cella@contig_tbl) %>% filter(clinical!='Healthy') %>% 
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
  # filter(primary!='CD8 Naive') %>% 
  ggplot(aes(primary,value,fill=secondary))+geom_col(position='stack')+
  scale_x_discrete(breaks=names(table(h1$primary)),
                   labels=paste0(names(table(h1$primary)),'\n(',as.character(table(h1$primary)),')'))+
  xlab('Primary (# of cells)')+theme_minimal()
  # don't have scaling worked out yet, and this is overall
  # but treg lot of overlap cd4act cd8, cd8exh-act strong, cd4 act-exh?
  # lot more treg and cd4 act in nlr, low treg in l, 
table(h1) %>% t() %>% as.data.frame.matrix() %>% 
  mutate(across(.fns=proportions)) %>%
  pheatmap::pheatmap(
    cluster_rows=F,cluster_cols=F,
    breaks=c(seq(0,0.4,length.out=100),0.6),
    color=colorRampPalette(RColorBrewer::brewer.pal(7,'Purples'))(100))

# lazy number of types per clonotype by clinical
cella@contig_tbl %>% filter(chain=='TRB') %>% 
  group_by(patient,timepoint,cdr3) %>%
  summarize(ntypes=if_else(length(na.omit(cluster_coarse))==0,NA_integer_,
                           length(unique(na.omit(cluster_coarse))))) %>% 
  na.omit() %>% add_clinical() %>% 
  mutate(ntypes=factor(ntypes)) %>% 
  group_by(patient,timepoint,clinical,ntypes) %>% 
  summarize(count=length(ntypes)) %>% ungroup() %>% 
  group_by(patient,timepoint,clinical) %>% 
  summarize(prop=proportions(count),ntypes=ntypes) %>% ungroup() %>% 
  filter(ntypes!='1') %>% mutate(ntypes=if_else(ntypes=='2','2','>2')) %>% #this line to just compare 2 and >2 number of matching types
  ggplot(aes(ntypes,prop))+geom_boxplot(aes(fill=clinical))+scale_y_sqrt()
  # not much clinical diff, maybe a bit more single in nr
  # more 2+ in nr, more 2 in nlr/maybe l; the high nr outliers are 17


# TM/BM diversity (superseded) ----------------------------------------------------------

# TM pbmc clonality
repDiversity(map(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
  x[!is.na(x$expTil),] #TM pbmc
}),.method='inv.simp') %>% 
  # as.data.frame() %>% select(Value=V1) %>% rownames_to_column('Sample') %>%
  separate(Sample,c('patient','timepoint')) %>% add_clinical() %>% 
  # ggplot(aes(reorder(interaction(patient,timepoint),desc(Value)),
  #            Value,fill=clinical))+
  # geom_col()+theme(axis.text.x=element_text(angle=90))
  plotBoxplot()+scale_y_log10()
  # ggplot(aes(timepoint,Value,fill=clinical))+geom_boxplot()
  # hard to really say any diff here
# regressing on nclonotypes - this is wrong, should regress on nseqs and then there's no difference - but why is it so strong for clonotypes, what does that mean
h1 <- map_dbl(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
  # x[!is.na(x$expTil),] %>% .$Clones %>% sum()
  nrow(x[!is.na(x$expTil),])
})
h1 <- repDiversity(
  map(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
  x[!is.na(x$expTil),] #TM pbmc
}),.method='inv.simp') %>% 
  # as.data.frame() %>% select(Value=Value) %>% rownames_to_column('Sample') %>%
  mutate(count=h1) %>% filter(count>1) %>% #remove 23-1
  separate(Sample,c('patient','timepoint')) %>% add_clinical()
summary(lm(Value~count,h1)) #p=2e-6, adj r^2=0.35
ggplot(h1,aes(count,Value))+geom_point(aes(color=clinical))+
  scale_y_log10()+scale_x_log10()+geom_smooth(method='lm')
mutate(h1,resid=lm(Value~count,h1)$residuals) %>% 
  # filter(timepoint=='5') %>%  #seems pretty valid at all timepoints
  ggplot(aes(reorder(interaction(patient,timepoint),desc(resid)),
             resid,fill=clinical))+
  geom_col()+theme(axis.text.x=element_text(angle=90))
  # holds for all timepoints
  #### what if we compared this to overall div of each
  # test diff div metrics for least signif corr between count and value
    # gini.simp 8.3F, inv.simp 29F, chao1 267F, gini 3.9F, d50 34F, hill1 92F, hill2 29F, hill3 16F, hill4 13F, hill5 11F
    # with resid: inv.simp-,gini+,gini.simp-,chao1x,d50-,all hill-

# BM TIL diversity
  #### somehow scale by the diversity of overall?
repDiversity(map(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
  x[!is.na(x$expTil),] %>% mutate(Clones=expTil)
}),.method='inv.simp') %>% 
  separate(Sample,c('patient','timepoint')) %>% add_clinical() %>% 
  plotBoxplot() # no real diff, 3/5 are l>nlr>nr
  # ggplot(aes(reorder(patient,desc(Value)),Value,fill=clinical))+geom_col()
h1 <- map_dbl(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
  # x[!is.na(x$expTil),] %>% .$Clones %>% sum()
  nrow(x[!is.na(x$expTil),])
})
h1 <- repDiversity(
  map(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
    x[!is.na(x$expTil),] #TM pbmc
  }),.method='inv.simp') %>% 
  # as.data.frame() %>% select(Value=V1) %>% rownames_to_column('Sample') %>%
  mutate(count=h1) %>% filter(count>1) %>% #remove 23-1
  separate(Sample,c('patient','timepoint')) %>% add_clinical()
summary(lm(Value~count,h1)) #p=2e-6, adj r^2=0.35
ggplot(h1,aes(count,Value))+geom_point(aes(color=clinical))+
  scale_y_log10()+scale_x_log10()+geom_smooth(method='lm')
mutate(h1,resid=lm(Value~count,h1)$residuals) %>% 
  ggplot(aes(reorder(interaction(patient,timepoint),desc(resid)),
             resid,fill=clinical))+
  geom_col()+theme(axis.text.x=element_text(angle=90))
  # same kind of patterns as TM pbmc - inconclusive overall, negative for most resid, positive on gini

# aight fuck it lazy 821 plot time
map_dfr(names(rep4$data[rep4$meta$timepoint %in% c('1')]),function(x){
  y <- na.omit(rep4$data[[x]]$expTil)
  p <- str_extract(x,'.*(?=-)')
  tp <- str_extract(x,'.$')
  data.frame('patient'=p,'timepoint'=tp,
             'high'=mean(y>=8),'low'=mean((y>1)&(y<8)),'singlet'=mean(y==1))
}) %>% add_clinical() %>% 
  pivot_longer(c(high,low,singlet),names_to='expansion') %>% 
  plotBoxplot('value')+facet_wrap(vars(expansion))
  # yeah still shows a lot more high in lynch, but not as pretty
  #### do like a gradient of color in a column but idk how

# 821 plot for BM expansion of BM tils
map_dfr(names(rep4$data[rep4$meta$timepoint %in% c('1')]),function(x){
  y <- rep4$data[[x]]$Clones[!is.na(rep4$data[[x]]$expTil)]
  p <- str_extract(x,'.*(?=-)')
  tp <- str_extract(x,'.$')
  data.frame('patient'=p,'timepoint'=tp,
             'high'=mean(y>=8),'low'=mean((y>1)&(y<8)),'singlet'=mean(y==1))
}) %>% add_clinical() %>% 
  pivot_longer(c(high,low,singlet),names_to='expansion') %>% 
  plotBoxplot('value')+facet_wrap(vars(expansion))
  # seems like actually more clonal nr, more singlet in responder

# does the earlier regression stuff indicate lower evenness in responders?
h1 <- map_dfr(rep4$data[(rep4$meta$timepoint %in% c('1'))&
                          (!(names(rep4$data) %in% c('23-1','20-5')))],function(x){
  x <- x[!is.na(x$expTil),]
  if(any(duplicated(x$CDR3.aa))) warning('Duplicated CDR3aa not accounted for')
  data.frame('div'=vegan::diversity(x$Clones,'shannon'),
             'nClones'=nrow(x),'nSeqs'=sum(x$Clones),
             'patient'=x$patient[1],'timepoint'=x$timepoint[1])
}) %>% add_clinical()
plotBoxplot(h1,'div')
summary(lm(div~nClones,h1))
ggplot(h1,aes(nClones,div))+geom_point(aes(color=clinical))+
  scale_y_log10()+scale_x_log10()+geom_smooth(method='lm')
mutate(h1,resid=lm(div~nClones,h1)$residuals) %>% 
  ggplot(aes(reorder(interaction(patient,timepoint),desc(resid)),
             resid,fill=clinical))+
  geom_col()+theme(axis.text.x=element_text(angle=90))
  # no because there isn't even a regression difference in diversity


# predictive diversity -----------------------------------------------------

# tp 0/1 diversity
rep4$data[rep4$meta$timepoint=='1'] %>% 
  # map(~filter(.x,cluster_fine=='CD8 TEMRA')) %>%
  repDiversity(.method='inv.simp') %>%
  mutate(patient=str_extract(Sample,'.*(?=-)')) %>% 
  add_clinical() %>% plotBoxplot()+scale_y_log10()+
  # ggrepel::geom_text_repel(aes(label=patient),max.overlaps=25)+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1))
  # lower div lynch for 1, holds for residualized version
  # l>nlr>nr div for 0


# predictive cell counts ----------------------------------------------------

# cell counts
  # for ldat: change dat in load and do cluster_fine and then set to res.0.9
h1 <- (dat$timepoint=='1')|(dat$clinical=='Healthy')
table(dat$cluster_coarse[h1],dat$clinical[h1],dat$patient[h1]) %>% 
  as.data.frame() %>% 
  set_names(c('type','clinical','patient','freq')) %>% 
  filter(freq>0) %>% group_split(patient,clinical) %>% 
  map_dfr(~mutate(.x,prop=proportions(freq))) %>% 
  plotBoxplot('prop')+facet_wrap(vars(type),scales='free')+
  theme(axis.text.x=element_blank())

# cell counts cluster
h1 <- (dat$timepoint=='1')#|(dat$clinical=='Healthy')
h1 <- table(dat$integrated_snn_res.1.3[h1],dat$clinical[h1],dat$patient[h1]) %>% 
  as.data.frame() %>% 
  set_names(c('type','clinical','patient','freq')) %>% 
  filter(freq>0) %>% group_split(patient,clinical) %>% 
  map_dfr(~mutate(.x,prop=proportions(freq)))
plotBoxplot(h1,'prop')+facet_wrap(vars(type),scales='free')+
  theme(axis.text.x=element_blank())
  # got some interesting ones, 25 is just 14, 15 strong 23/5, 23-1 is just 6 cells, 21-1 is 200 cells
# specific clusts - 3 (l high - temra right),8 (nlr lo - cd8 naive), 12 (nlr hi, cd4 activated ish, center top), 16 (nlr hi - dead?), 18(nlr hi,l lo, cd4 naive),21 (nlr hi - cd4 early?), 24 (l lo,cd8-4 transition), 29 (l lo - left gdt), 
# ldat clusts - 4(nlr lo-cd14mono),9(nlr hi-cd16 mono),19(l lo-cdc),21(nlr hi-cd16nk.2),22(l lo-cd56 nk),26(nlr hi-pdc),28(l lo,nlr hi-mono??),29(l hi-plasma),31(l lo-progenitor),32(l lo-b/mono)
filter(h1,type=='29') %>% plotBoxplot('prop')#+ggrepel::geom_text_repel(aes(label=patient))


# temporal TIL overlap --------------------------------------------------------

with15 <- group_by(rep4$meta,patient) %>% 
  summarize(a=if_else(all(c('1','5') %in% timepoint),patient[1],NA_character_)) %>%
  na.omit() %>% .$patient
map_dbl(with15,function(p){
  repOverlap(list(rep4$data[[paste0(p,'-1')]],rep4$data[[paste0(p,'-5')]]),.method='jaccard',.col='aa')
}) %>% set_names(with15) %>% 
  as.data.frame() %>% rownames_to_column('patient') %>% add_clinical() %>% 
  ggplot(aes(reorder(patient,desc(`.`)),`.`,fill=clinical))+geom_col()
# yo jaccard ignores proportions/clones, just based off cdr3s
map_dbl(with15,function(p){# this for TM
  list(
    rep4$data[[paste0(p,'-1')]] %>% filter(!is.na(expTil)),
    rep4$data[[paste0(p,'-5')]] %>% filter(!is.na(expTil))
  ) %>% repOverlap(.method='jaccard',.col='aa')
}) %>% set_names(with15) %>% 
  as.data.frame() %>% rownames_to_column('patient') %>% add_clinical() %>% 
  ggplot(aes(reorder(patient,desc(`.`)),`.`,fill=clinical))+geom_col()
# change all 5s to 3s for timepoint 3 stuff


# cd8s count/clonality -----------------------------------------------------------

# proportions of t cells
filter(dat,cluster_coarse!='Junk') %>% select(cluster_coarse,clinical) %>% 
  table() %>% as.data.frame.matrix() %>% rownames_to_column('cluster_coarse') %>% 
  mutate(across(-cluster_coarse,.fns=proportions)) %>% 
  filter(cluster_coarse=='CD8 Activated') %>% 
  pivot_longer(-cluster_coarse) %>% 
  ggplot(aes(reorder(name,value),value,fill=name))+geom_col()+theme_void()
# proportion of t cells split by fine type
filter(dat,cluster_coarse!='Junk') %>% select(cluster_fine,clinical) %>% 
  table() %>% as.data.frame.matrix() %>% rownames_to_column('clust') %>% 
  mutate(across(-clust,.fns=proportions)) %>% 
  filter(clust %in% c('CD8 Effector Memory','CD8 TEMRA')) %>% 
  pivot_longer(-clust) %>% 
  ggplot(aes(name,value,fill=factor(clust,c('CD8 TEMRA','CD8 Effector Memory'))))+geom_col(position='stack')+theme(legend.title=element_blank())
# proportions of overall
h1 <- readRDS('data_sc/meta_ldat.rds')
h1 <- table(h1$clinical) %>% as.vector()
h2 <- select(dat,cluster_coarse,clinical) %>% table()
t(t(h2)/as.vector(h1)) %>% as.data.frame() %>% 
  filter(cluster_coarse=='CD8 Activated') %>% 
  ggplot(aes(reorder(clinical,Freq),Freq,fill=clinical))+geom_col()+theme_void()
# proportions of t cells split by patient
group_split(dat,patient,timepoint) %>% 
  map(~table(.x$cluster_coarse) %>% as.data.frame() %>% 
        set_names('type',paste0(.x$patient[1],'-',.x$timepoint[1]))) %>% 
  reduce(left_join,by='type') %>% 
  mutate(across(-type,~proportions(replace_na(.x,0)))) %>% 
  pivot_longer(-type,values_to='prop') %>% 
  separate(name,c('patient','timepoint')) %>% 
  filter(!(interaction(patient,timepoint) %in% c('20.5','23.1')),
         type=='CD8 Activated') %>% 
  add_clinical() %>% plotBoxplot('prop')+ggrepel::geom_text_repel(aes(label=interaction(patient,timepoint)),max.overlaps=100)
# proportion of overall split by patient
h1 <- readRDS('data_sc/meta_ldat.rds')
h1 <- table(h1$patient,h1$timepoint) %>% as.data.frame() %>% 
  filter(Freq!=0) %>% arrange(as.character(Var1),Var2) %>% .$Freq
h2 <- group_split(dat,patient,timepoint) %>% 
  map(~table(.x$cluster_coarse) %>% as.data.frame() %>% 
        set_names('type',paste0(.x$patient[1],'-',.x$timepoint[1]))) %>% 
  reduce(left_join,by='type')
t(t(h2[,-1])/h1) %>% as.data.frame() %>% mutate(type=h2$type) %>% 
  pivot_longer(-type,values_to='prop') %>% 
  separate(name,c('patient','timepoint')) %>% 
  filter(!(interaction(patient,timepoint) %in% c('20.5','23.1')),
         type=='CD8 Activated') %>% 
  add_clinical() %>% plotBoxplot('prop')+ggrepel::geom_text_repel(aes(label=interaction(patient,timepoint)),max.overlaps=100)
  

# clonality stuff
# overall pretreatment
rep4$data %>% 
  repDiversity(.method='inv.simp') %>%
  mutate(patient=str_extract(Sample,'.*(?=-)'),
         timepoint=str_extract(Sample,'.$')) %>% 
  add_clinical() %>% mutate(clinical=factor(clinical,c('Lynch-like responder','Nonlynch-like responder','Nonresponder','Healthy'))) %>% 
  filter(timepoint %in% c(1,'h')) %>%
  ggplot(aes(clinical,Value,fill=clinical,color=clinical))+
  geom_boxplot(lwd=1.5)+scale_y_log10()+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1))+
  theme_prism(axis_text_angle=45,base_size=12)+
  scale_fill_manual(values=cols_clinical)+
  scale_color_manual(values=cols_clinical_dark)+
  ylab('Inverse simpson index')
# posttreatment
rep4$data %>% 
  # map(~filter(.x,cluster_fine=='CD8 TEMRA')) %>%
  repDiversity(.method='inv.simp') %>%
  mutate(patient=str_extract(Sample,'.*(?=-)'),
         timepoint=str_extract(Sample,'.$')) %>% 
  add_clinical() %>% mutate(clinical=factor(clinical,c('Lynch-like responder','Nonlynch-like responder','Nonresponder','Healthy'))) %>% 
  filter(timepoint %in% c('3','5','h')) %>% 
  ggplot(aes(clinical,Value,fill=clinical,color=clinical))+
  geom_boxplot(lwd=1.5)+scale_y_log10()+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1))+
  theme_prism(axis_text_angle=45,base_size=12)+
  scale_fill_manual(values=cols_clinical)+
  scale_color_manual(values=cols_clinical_dark)+
  ylab('Inverse simpson index')
# checking that isn't just seq number question
h1 <- repDiversity(rep4$data,.method='inv.simp') %>%
  mutate(patient=str_extract(Sample,'.*(?=-)'),
         timepoint=str_extract(Sample,'.$')) %>% add_clinical()
h1$seqs <- map_dbl(rep4$data,~sum(.x$Clones))
filter(h1,clinical=='Healthy') %>% ggplot(aes(seqs,Value))+geom_point(aes(color=clinical,pch=timepoint))+geom_smooth(method='lm')+scale_y_log10()+scale_x_log10()
filter(h1,timepoint %in% c('1')) %>% ggplot(aes(seqs,Value))+geom_point(aes(color=clinical,pch=timepoint))+geom_smooth(method='lm')+scale_y_log10()+scale_x_log10()
filter(h1,timepoint %in% c('3','5')) %>% ggplot(aes(seqs,Value))+geom_point(aes(color=clinical,pch=timepoint))+geom_smooth(method='lm')+scale_y_log10()+scale_x_log10()
filter(h1,timepoint %in% c('0','7')) %>% ggplot(aes(seqs,Value))+geom_point(aes(color=clinical,pch=timepoint))+geom_smooth(method='lm')+scale_y_log10()+scale_x_log10()
# temra specific
rep4$data[(rep4$meta$timepoint %in% c(1,3,5))&
            (!names(rep4$data) %in% c('20-5','23-1'))] %>% 
  map(~filter(.x,cluster_fine=='CD8 TEMRA')) %>%
  repDiversity(.method='inv.simp') %>%
  mutate(patient=str_extract(Sample,'.*(?=-)'),
         timepoint=str_extract(Sample,'.$')) %>% 
  add_clinical() %>% mutate(clinical=factor(clinical,c('Lynch-like responder','Nonlynch-like responder','Nonresponder','Healthy'))) %>% 
  ggplot(aes(clinical,Value,fill=clinical,color=clinical))+
  geom_boxplot(lwd=1.5)+scale_y_log10()+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1))+
  theme_prism(axis_text_angle=45,base_size=12)+
  scale_fill_manual(values=cols_clinical)+
  scale_color_manual(values=cols_clinical_dark)+
  ylab('Inverse simpson index')+facet_wrap(~timepoint)


# TEMRA distinctness -------------------------------------------------------

#' want to do an overlap among cell types, averaged over patients, as circos
#' too complex, try overlap means of everything to temra, radar plot
names(rep4$data)[(rep4$meta$timepoint %in% c(1,3,5,'h'))&
    (!names(rep4$data) %in% c('20-5','23-1'))] %>% 
  map(function(n){
    #' loop over data and split by cell type
      #' needs to be able to fail
    #' loop over cell types and run jaccard
    h1 <- rep4$data[[n]] %>% filter(!is.na(cluster_fine)) %>% 
      mutate(temp=if_else(cluster_fine=='CD8 TEMRA','CD8 TEMRA','Other')) %>% 
      group_split(temp)
    names(h1) <- map_chr(h1,~.x$temp[1])
    # if('CD8 TEMRA' %in% names(h1)){
    #   map_dfr(names(h1)[names(h1)!='CD8 TEMRA'],function(m){
    #     print(names(h1))
    #     h2 <- repOverlap(list(h1[[m]],h1[['CD8 TEMRA']]),'jaccard') %>% 
    #       as.numeric()
    #     data.frame('patient'=str_extract(n,'.*(?=-)'),
    #                'timepoint'=str_extract(n,'.$'),
    #                'type'=m,'value'=h2)
    #   })
    # } else data.frame('patient'=str_extract(n,'.*(?=-)'),
    #                   'timepoint'=str_extract(n,'.$'),
    #                   'type'=NA,'value'=NA)
    map(h1,~nrow(.x))
  })
# temra v all split by tp/p, manual morisita
limitSize <- 5
filter(cella@cell_tbl,!is.na(cluster_fine)) %>% 
  group_split(patient,timepoint) %>% 
  map_dbl(function(x){
    if((sum(x$cluster_fine=='CD8 TEMRA')<limitSize)|(sum(x$cluster_fine!='CD8 TEMRA')<limitSize)) NA_real_
    else{
      h1 <- x$cdr3[x$cluster_fine=='CD8 TEMRA'] %>% 
        table() %>% as.data.frame() %>% 
        set_colnames(c('cdr3','temra'))
      h2 <- x$cdr3[x$cluster_fine!='CD8 TEMRA'] %>% 
        table() %>% as.data.frame() %>% 
        set_colnames(c('cdr3','other'))
      full_join(h1,h2,by='cdr3') %>% 
        replace_na(list('temra'=0,'other'=0)) %>% 
        select(-cdr3) %>% t() %>% 
        vegan::vegdist(method='morisita',binary=F)
        # vegan::vegdist(method='jaccard',binary=T)
    }
  }) %>% as.data.frame() %>% transmute(value=1-`.`) %>% 
  #output is dissimilarity (1 is no overlap), this converts output to similarity
  cbind(filter(cella@cell_tbl,!is.na(cluster_fine)) %>% 
          group_by(patient,timepoint) %>% group_keys()) %>% 
  add_clinical() %>% na.omit() %>% 
  plotBoxplot('value')
  # group_by(clinical) %>% 
  # summarize(avg=median(value),err=sd(value)/sqrt(length(value))) %>% 
  # ggplot(aes(clinical,avg,fill=clinical))+geom_col()+
  # geom_errorbar(aes(ymax=err+avg,ymin=avg-err),width=0.2,lwd=1)+
  # theme_prism()
# temra v all types
limitSize <- 5 #minimum n_cdr3 in each group to run overlap
fineTypes <- unique(cella@cell_tbl$cluster_fine)
fineTypes <- fineTypes[!is.na(fineTypes)&(fineTypes!='CD8 TEMRA')]
h1 <- filter(cella@cell_tbl,!is.na(cluster_fine),!(is.na(cdr3))) %>% 
  group_split(patient,timepoint) %>% 
  map_dfr(function(x){
    map_dbl(fineTypes,function(type2){
      if((sum(x$cluster_fine=='CD8 TEMRA')<limitSize)|
         (sum(x$cluster_fine==type2)<limitSize)) NA_real_
      else{
        h1 <- x$cdr3[x$cluster_fine=='CD8 TEMRA'] %>% 
          table() %>% as.data.frame() %>% 
          set_colnames(c('cdr3','temra'))
        h2 <- x$cdr3[x$cluster_fine==type2] %>% 
          table() %>% as.data.frame() %>% 
          set_colnames(c('cdr3','other'))
        full_join(h1,h2,by='cdr3') %>% 
          replace_na(list('temra'=0,'other'=0)) %>% 
          select(-cdr3) %>% t() %>% 
          vegan::vegdist(method='jaccard',binary=T)
      }
    }) %>% matrix(nrow=1) %>% subtract(e1=1,e2=.) %>% 
      as.data.frame() %>% set_names(fineTypes) %>% 
      mutate(patient=x$patient[1],timepoint=x$timepoint[1])
  }) %>% add_clinical()
plotBoxplot(h1,'CD8 Effector Memory')
pivot_longer(h1,-c(patient,clinical,timepoint),
             names_to='type',values_to='similarity') %>% 
  ggplot(aes(y=similarity,fill=clinical))+geom_boxplot()+
  facet_wrap(~type,scales='free')
pivot_longer(h1,-c(patient,clinical,timepoint),
             names_to='type',values_to='similarity') %>% 
  filter(type %in% c('CD8 Effector Memory','Cycling','Th1','Th2','Treg')) %>% 
  group_by(clinical,type) %>%
  summarize(avg=mean(similarity,na.rm=T),
            err=sd(similarity,na.rm=T)/sqrt(length(na.omit(similarity)))) %>%
  ggplot(aes(clinical,avg,fill=clinical,color=clinical))+
  geom_col(position=position_dodge(width=1))+
  geom_errorbar(aes(ymax=err+avg,ymin=avg-err),
                width=0.2,lwd=0.8,position=position_dodge(width=1))+
  facet_wrap(vars(type),scales='free')+
  ylab('Jaccard similarity')+theme_prism()+
  theme(axis.text.x=element_blank())+
  scale_fill_manual(values=cols_clinical)+
  scale_color_manual(values=cols_clinical_dark)
  # still not entirely sure about the mechanics of vegdist
  # could still be a size effect of the temra, haven't done downsampling
  # also mean or median for the column version?


# change in diversity -----------------------------------------------------

# change in diversity
h1 <- rep4$data[(!(names(rep4$data) %in% c('20-5','23-1')))&
            (rep4$meta$timepoint %in% c(1,3,5))] %>% 
  repDiversity('inv.simp') %>% mutate(Value=log(Value)) %>% 
  separate(Sample,c('patient','timepoint')) %>% 
  pivot_wider(id_cols=patient,names_from=timepoint,values_from=Value)
  # 1-3-5 line plot
na.omit(h1) %>% mutate(`3`=(`3`-`1`)/`1`,`5`=(`5`-`1`)/`1`,`1`=0) %>% 
  pivot_longer(-patient) %>% add_clinical() %>% ggplot(aes(name,value))+
  geom_line(aes(group=patient,color=clinical))
  # rainfall of div change
mutate(h1,`3-1`=(`3`-`1`)/`1`,`5-1`=(`5`-`1`)/`1`,`5-3`=(`5`-`3`)/`3`) %>% 
  add_clinical() %>% 
  ggplot(aes(reorder(patient,desc(`5-1`)),`5-1`,fill=clinical))+geom_col()+scale_fill_manual(values=cols_clinical)+theme_prism()+theme(axis.text.x=element_text(angle=90))+labs(y='Change in diversity')
# change in TM diversity
h1 <- rep4$data[(!(names(rep4$data) %in% c('20-5','23-1')))&
            (rep4$meta$timepoint %in% c(1,3,5))]
h1 <- map_dbl(h1,~repDiversity(mutate(.x,Clones=expTil),'inv.simp') %>% 
                as.numeric()) %>% 
  data.frame('Value'=.,'Sample'=names(h1)) %>% 
  mutate(Value=log(Value)) %>% separate(Sample,c('patient','timepoint')) %>% 
  pivot_wider(id_cols=patient,names_from=timepoint,values_from=Value)
  # 1-3-5 line plot
na.omit(h1) %>% mutate(`3`=(`3`-`1`)/`1`,`5`=(`5`-`1`)/`1`,`1`=0) %>% 
  pivot_longer(-patient) %>% add_clinical() %>% ggplot(aes(name,value))+
  geom_line(aes(group=patient,color=clinical))
  # rainfall of div change
mutate(h1,`3-1`=(`3`-`1`)/`1`,`5-1`=(`5`-`1`)/`1`,`5-3`=(`5`-`3`)/`3`) %>% 
  add_clinical() %>% 
  ggplot(aes(reorder(patient,desc(`5-1`)),`5-1`,fill=clinical))+geom_col()+scale_fill_manual(values=cols_clinical)+theme_prism()+theme(axis.text.x=element_text(angle=90))+labs(y='Change in diversity')



# proportion of clonotypes over size n ------------------------------------

# trying to determine what prop of clonotypes (or seqs?) are in clonotypes over size n
clone_prop <- function(df,n,margin='seqs'){
  if(margin=='seqs') sum(filter(df,Clones>=n)$Proportion)
  else if(margin=='clones') mean(df$Clones>=n)
  else stop('invalid margin')
}
# now run some like exponential number series for n for each to generate a curve of prop of clonotypes/seqs vs n
clone_curves <- function(dfs,margin='seqs',max1=20){
  clone_curve <- function(df,margin='seqs',max1=20){
    map_dfr(1:max1,function(i){
      data.frame('n'=i,'prop'=clone_prop(df=df,n=i,margin=margin))
    })
  }
  if(any(class(dfs)=='data.frame')){
    out <- clone_curve(df=dfs,margin=margin,max1=max1) %>% 
      mutate(patient='?',clinical='?',timepoint='?')
  }
  else if(any(class(dfs)=='list')){
    out <- map_dfr(dfs,.id='name',~clone_curve(df=.x,margin=margin,max1=max1)) %>% 
      separate(name,c('patient','timepoint')) %>% add_clinical()
  }
  label_df <- filter(out,n==2) %>% select(patient,timepoint,clinical,prop)
  ggplot(out,aes(n,prop,color=clinical))+
    geom_line(aes(group=interaction(patient,timepoint)))+
    scale_x_log10()+ylim(c(0,1))+
    labs(y=paste('Proportion of',margin),x='Size threshold')+
    ggrepel::geom_label_repel(data=label_df,aes(2,prop,color=clinical,
                                           label=interaction(patient,timepoint)))
}
clone_curves(rep4$data,'seqs',20)

clone_curve2 <- function(dfs,margin='seqs',max1=20){
  clone_curve <- function(df,margin='seqs',max1=20){
    map_dfr(1:max1,function(i){
      data.frame('n'=i,'prop'=clone_prop(df=df,n=i,margin=margin))
    })
  }
  if(any(class(dfs)=='data.frame')){
    out <- clone_curve(df=dfs,margin=margin,max1=max1) %>% 
      mutate(patient='?',clinical='?',timepoint='?')
  }
  else if(any(class(dfs)=='list')){
    out <- map_dfr(dfs,.id='name',~clone_curve(df=.x,margin=margin,max1=max1)) %>% 
      separate(name,c('patient','timepoint')) %>% add_clinical()
  }
  # we need to uhh
  add_count(out,clinical,n,name='count') %>% 
    group_by(clinical,n) %>%
    summarize(serr=sd(prop)/sqrt(count[1]),prop=mean(prop)) %>% ungroup() %>%
    ggplot(aes(n,prop,color=clinical))+
    geom_point(size=1.25)+geom_line(aes(group=clinical))+
    geom_errorbar(aes(ymin=prop-serr,ymax=prop+serr),width=0.2)+
    ylim(c(0,1))+
    labs(y=paste('Proportion of',margin),x='Size threshold')
}
clone_curve2(rep4$data[(rep4$meta$timepoint %in% c('1','3','5','h'))&
                         (!(names(rep4$data) %in% c('20-5','23-1')))])+
  ggprism::theme_prism()+scale_color_manual(values=cols_clinical)


# just diversity ----------------------------------------------------------

### overall diversity
h1 <- rep4$data[!(names(rep4$data) %in% c('20-5','23-1'))] %>% 
  repDiversity('inv.simp') %>% 
  separate(Sample,c('patient','timepoint')) %>% add_clinical()

## raw diversity
# pbmc overall
filter(h1,timepoint %in% c('1','3','5','h')) %>% 
  plotBoxplot() %>% prism()+scale_y_log10()
# pbmc pretreatment
filter(h1,timepoint %in% c('1','h')) %>% plotBoxplot() %>% prism()+scale_y_log10()
# pbmc posttreatment
filter(h1,timepoint %in% c('3','5','h')) %>% 
  plotBoxplot() %>% prism()+scale_y_log10()
# pbmc timepoint 5
filter(h1,timepoint %in% c('5','h')) %>% plotBoxplot() %>% prism()+scale_y_log10()
# pbmc split timepoint box
filter(h1,timepoint %in% c('1','3','5')) %>% 
  plotBoxplot() %>% prism()+scale_y_log10()+facet_wrap(~timepoint)+
  theme(axis.text.x=element_blank())
# pbmc split timepoint line+errorbars
filter(h1,timepoint %in% c('1','3','5','h')) %>% 
  mutate(Value=log10(Value)) %>% #weird geom mean geom sterr stuff here
  group_by(clinical,timepoint) %>%
  summarize(serr=sd(Value)/sqrt(n()),Value=mean(Value)) %>% ungroup() %>%
  ggplot(aes(timepoint,10**Value,color=clinical))+
  geom_point(size=2.5,position=position_dodge(width=0.15))+
  geom_line(aes(group=clinical),position=position_dodge(width=0.15),lwd=1)+
  geom_errorbar(aes(ymin=10**(Value-serr),ymax=10**(Value+serr)),width=0.1,
                position=position_dodge(width=0.15),lwd=1)+
  scale_y_log10()+theme_prism()+scale_color_manual(values=cols_clinical)
# pbmc split timepoint line+errorbars (3-5 combined)
filter(h1,timepoint %in% c('1','3','5')) %>% 
  mutate(Value=log10(Value),timepoint=if_else(timepoint=='1','pre','post') %>% 
           factor(levels=c('pre','post'))) %>%
  group_by(clinical,timepoint) %>%
  summarize(serr=sd(Value)/sqrt(n()),Value=mean(Value)) %>% ungroup() %>%
  ggplot(aes(timepoint,10**Value,color=clinical))+
  geom_point(size=2.5,position=position_dodge(width=0.15))+
  geom_line(aes(group=clinical),position=position_dodge(width=0.15),lwd=1)+
  geom_errorbar(aes(ymin=10**(Value-serr),ymax=10**(Value+serr)),width=0.1,
                position=position_dodge(width=0.15),lwd=1)+
  scale_y_log10()+theme_prism()+scale_color_manual(values=cols_clinical)
# TIL
filter(h1,timepoint=='0') %>% plotBoxplot() %>% prism()+scale_y_log10()

## Change in diversity
# pbmc
h2 <- filter(h1,timepoint %in% c('1','3','5')) %>% 
  mutate(Value=log10(Value)) %>% #Values already log scaled
  pivot_wider(id_cols=patient,names_from=timepoint,values_from=Value)
# 1-3-5 line plot (proportional)
na.omit(h2) %>% mutate(`3`=(`3`-`1`)/`1`,`5`=(`5`-`1`)/`1`,`1`=0) %>% 
  pivot_longer(-patient) %>% add_clinical() %>% ggplot(aes(name,value))+
  geom_line(aes(group=patient,color=clinical),lwd=1.2)+
  labs(y='Proportional change',x='Timepoint',title='Diversity over time')+
  theme_prism()+scale_color_manual(values=cols_clinical)
# rainfall of div change (proportional)
mutate(h2,`3-1`=(`3`-`1`)/`1`,`5-1`=(`5`-`1`)/`1`,`5-3`=(`5`-`3`)/`3`) %>% 
  add_clinical() %>% filter(!is.na(`5-1`)) %>% 
  ggplot(aes(reorder(patient,desc(`5-1`)),`5-1`,fill=clinical)) %>% 
  prism()+geom_col()+theme(axis.text.x=element_text(angle=90))+
  labs(y='Change in diversity',x='Patient')
# errobars of div change (proportional)
mutate(h2,`3-1`=(`3`-`1`)/`1`,`5-1`=(`5`-`1`)/`1`,`5-3`=(`5`-`3`)/`3`) %>% 
  add_clinical() %>% filter(!is.na(`5-1`)) %>% 
  group_by(clinical) %>% summarize(Value=mean(`5-1`),serr=sd(`5-1`)/sqrt(n())) %>% 
  ggplot(aes(clinical,Value,fill=clinical,color=clinical)) %>% prism()+geom_col()+
  geom_errorbar(aes(ymin=Value-serr,ymax=Value+serr),width=0.2,lwd=1)
# errorbars of div change for each timepoint combination (proportional)
mutate(h2,`3-1`=(`3`-`1`)/`1`,`5-1`=(`5`-`1`)/`1`,`5-3`=(`5`-`3`)/`3`) %>% 
  select(-`1`,-`3`,-`5`) %>% pivot_longer(-patient) %>% add_clinical() %>% 
  group_by(clinical,name) %>% summarize(
    serr=sd(value,na.rm=T)/sqrt(sum(!is.na(value))),
    value=mean(value,na.rm=T)) %>% 
  ggplot(aes(clinical,value,color=clinical))+
  geom_point(aes(pch=name),position=position_dodge(width=0.3),size=2.5)+
  geom_errorbar(aes(ymin=value-serr,ymax=value+serr,group=name),width=0.2,lwd=1,
                position=position_dodge(width=0.3))+
  theme_prism()+scale_color_manual(values=cols_clinical)+
  theme(axis.text.x=element_blank())+geom_hline(yintercept=0)
# rainfall of div change (raw)
mutate(h2,`3-1`=(`3`-`1`),`5-1`=(`5`-`1`),`5-3`=(`5`-`3`)) %>% 
  add_clinical() %>% filter(!is.na(`5-1`)) %>% 
  ggplot(aes(reorder(patient,desc(`5-1`)),`5-1`,fill=clinical)) %>% 
  prism()+geom_col()+theme(axis.text.x=element_text(angle=90))+
  labs(y='Change in diversity',x='Patient')


### TM diversity
h1 <- rep4$data[(rep4$meta$timepoint %in% c(1,3,5))] %>% 
  map(~mutate(.x,Clones=expTil) %>% filter(!is.na(Clones)))
h2 <- map_dbl(h1,nrow) %T>% print()
h1 <- h1[h2>=10] #can set min clones cutoff
h1 %<>% repDiversity('inv.simp') %>% mutate(nclones=h2[Sample]) %>% 
  separate(Sample,c('patient','timepoint')) %>% add_clinical()

## raw diversity
# pbmc overall
filter(h1,timepoint %in% c('1','3','5','h')) %>% 
  plotBoxplot() %>% prism()+scale_y_log10()
# pbmc pretreatment
filter(h1,timepoint %in% c('1','h')) %>% plotBoxplot() %>% prism()+scale_y_log10()
# pbmc posttreatment
filter(h1,timepoint %in% c('3','5','h')) %>% 
  plotBoxplot() %>% prism()+scale_y_log10()
# pbmc timepoint 5
filter(h1,timepoint %in% c('5','h')) %>% plotBoxplot() %>% prism()+scale_y_log10()
# pbmc split timepoint box
filter(h1,timepoint %in% c('1','3','5')) %>% 
  plotBoxplot() %>% prism()+scale_y_log10()+facet_wrap(~timepoint)+
  theme(axis.text.x=element_blank())
# pbmc split timepoint line+errorbars
filter(h1,timepoint %in% c('1','3','5')) %>% 
  mutate(Value=log10(Value)) %>% #weird geom mean geom sterr stuff here
  group_by(clinical,timepoint) %>%
  summarize(serr=sd(Value)/sqrt(n()),Value=mean(Value)) %>% ungroup() %>%
  ggplot(aes(timepoint,10**Value,color=clinical))+
  geom_point(size=2.5,position=position_dodge(width=0.15))+
  geom_line(aes(group=clinical),position=position_dodge(width=0.15),lwd=1)+
  geom_errorbar(aes(ymin=10**(Value-serr),ymax=10**(Value+serr)),width=0.1,
                position=position_dodge(width=0.15),lwd=1)+
  scale_y_log10()+theme_prism()+scale_color_manual(values=cols_clinical)
# pbmc split timepoint line+errorbars (3-5 combined)
filter(h1,timepoint %in% c('1','3','5')) %>% 
  mutate(Value=log10(Value),timepoint=if_else(timepoint=='1','pre','post') %>% 
           factor(levels=c('pre','post'))) %>%
  group_by(clinical,timepoint) %>%
  summarize(serr=sd(Value)/sqrt(n()),Value=mean(Value)) %>% ungroup() %>%
  ggplot(aes(timepoint,10**Value,color=clinical))+
  geom_point(size=2.5,position=position_dodge(width=0.15))+
  geom_line(aes(group=clinical),position=position_dodge(width=0.15),lwd=1)+
  geom_errorbar(aes(ymin=10**(Value-serr),ymax=10**(Value+serr)),width=0.1,
                position=position_dodge(width=0.15),lwd=1)+
  scale_y_log10()+theme_prism()+scale_color_manual(values=cols_clinical)

## Change in diversity
# pbmc
h1 <- filter(h1,timepoint %in% c('1','3','5')) %>% 
  mutate(Value=log10(Value)) %>% #Values already log scaled
  pivot_wider(id_cols=patient,names_from=timepoint,values_from=Value)
# 1-3-5 line plot (proportional)
na.omit(h1) %>% mutate(`3`=(`3`-`1`)/`1`,`5`=(`5`-`1`)/`1`,`1`=0) %>% 
  pivot_longer(-patient) %>% add_clinical() %>% ggplot(aes(name,value))+
  geom_line(aes(group=patient,color=clinical),lwd=1.2)+
  labs(y='Proportional change',x='Timepoint',title='Diversity over time')+
  theme_prism()+scale_color_manual(values=cols_clinical)
# rainfall of div change (proportional)
mutate(h1,`3-1`=(`3`-`1`)/`1`,`5-1`=(`5`-`1`)/`1`,`5-3`=(`5`-`3`)/`3`) %>% 
  add_clinical() %>% filter(!is.na(`5-1`)) %>% 
  ggplot(aes(reorder(patient,desc(`5-1`)),`5-1`,fill=clinical)) %>% 
  prism()+geom_col()+theme(axis.text.x=element_text(angle=90))+
  labs(y='Change in diversity',x='Patient')
# errobars of div change (proportional)
mutate(h1,`3-1`=(`3`-`1`)/`1`,`5-1`=(`5`-`1`)/`1`,`5-3`=(`5`-`3`)/`3`) %>% 
  add_clinical() %>% filter(!is.na(`5-1`)) %>% 
  group_by(clinical) %>% summarize(Value=mean(`5-1`),serr=sd(`5-1`)/sqrt(n())) %>% 
  ggplot(aes(clinical,Value,fill=clinical,color=clinical)) %>% prism()+geom_col()+
  geom_errorbar(aes(ymin=Value-serr,ymax=Value+serr),width=0.2,lwd=1)
# errorbars of div change for each timepoint combination (proportional)
mutate(h1,`3-1`=(`3`-`1`)/`1`,`5-1`=(`5`-`1`)/`1`,`5-3`=(`5`-`3`)/`3`) %>% 
  select(-`1`,-`3`,-`5`) %>% pivot_longer(-patient) %>% add_clinical() %>% 
  group_by(clinical,name) %>% summarize(
    serr=sd(value,na.rm=T)/sqrt(sum(!is.na(value))),
    value=mean(value,na.rm=T)) %>% 
  ggplot(aes(clinical,value,color=clinical))+
  geom_point(aes(pch=name),position=position_dodge(width=0.3),size=2.5)+
  geom_errorbar(aes(ymin=value-serr,ymax=value+serr,group=name),width=0.2,lwd=1,
                position=position_dodge(width=0.3))+
  theme_prism()+scale_color_manual(values=cols_clinical)+
  theme(axis.text.x=element_blank())+geom_hline(yintercept=0)
# rainfall of div change (raw)
mutate(h1,`3-1`=(`3`-`1`),`5-1`=(`5`-`1`),`5-3`=(`5`-`3`)) %>% 
  add_clinical() %>% filter(!is.na(`5-1`)) %>% 
  ggplot(aes(reorder(patient,desc(`5-1`)),`5-1`,fill=clinical)) %>% 
  prism()+geom_col()+theme(axis.text.x=element_text(angle=90))+
  labs(y='Change in diversity',x='Patient')





# TM pbmc clonality
repDiversity(map(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
  x[!is.na(x$expTil),] #TM pbmc
}),.method='inv.simp') %>% 
  # as.data.frame() %>% select(Value=V1) %>% rownames_to_column('Sample') %>%
  separate(Sample,c('patient','timepoint')) %>% add_clinical() %>% 
  # ggplot(aes(reorder(interaction(patient,timepoint),desc(Value)),
  #            Value,fill=clinical))+
  # geom_col()+theme(axis.text.x=element_text(angle=90))
  plotBoxplot()+scale_y_log10()
# ggplot(aes(timepoint,Value,fill=clinical))+geom_boxplot()
# hard to really say any diff here
# regressing on nclonotypes - this is wrong, should regress on nseqs and then there's no difference - but why is it so strong for clonotypes, what does that mean
h1 <- map_dbl(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
  # x[!is.na(x$expTil),] %>% .$Clones %>% sum()
  nrow(x[!is.na(x$expTil),])
})
h1 <- repDiversity(
  map(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
    x[!is.na(x$expTil),] #TM pbmc
  }),.method='inv.simp') %>% 
  # as.data.frame() %>% select(Value=Value) %>% rownames_to_column('Sample') %>%
  mutate(count=h1) %>% filter(count>1) %>% #remove 23-1
  separate(Sample,c('patient','timepoint')) %>% add_clinical()
summary(lm(Value~count,h1)) #p=2e-6, adj r^2=0.35
ggplot(h1,aes(count,Value))+geom_point(aes(color=clinical))+
  scale_y_log10()+scale_x_log10()+geom_smooth(method='lm')
mutate(h1,resid=lm(Value~count,h1)$residuals) %>% 
  # filter(timepoint=='5') %>%  #seems pretty valid at all timepoints
  ggplot(aes(reorder(interaction(patient,timepoint),desc(resid)),
             resid,fill=clinical))+
  geom_col()+theme(axis.text.x=element_text(angle=90))
# holds for all timepoints
#### what if we compared this to overall div of each
# test diff div metrics for least signif corr between count and value
# gini.simp 8.3F, inv.simp 29F, chao1 267F, gini 3.9F, d50 34F, hill1 92F, hill2 29F, hill3 16F, hill4 13F, hill5 11F
# with resid: inv.simp-,gini+,gini.simp-,chao1x,d50-,all hill-

# BM TIL diversity
#### somehow scale by the diversity of overall?
repDiversity(map(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
  x[!is.na(x$expTil),] %>% mutate(Clones=expTil)
}),.method='inv.simp') %>% 
  separate(Sample,c('patient','timepoint')) %>% add_clinical() %>% 
  plotBoxplot() # no real diff, 3/5 are l>nlr>nr
# ggplot(aes(reorder(patient,desc(Value)),Value,fill=clinical))+geom_col()
h1 <- map_dbl(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
  # x[!is.na(x$expTil),] %>% .$Clones %>% sum()
  nrow(x[!is.na(x$expTil),])
})
h1 <- repDiversity(
  map(rep4$data[rep4$meta$timepoint %in% c('1')],function(x){
    x[!is.na(x$expTil),] #TM pbmc
  }),.method='inv.simp') %>% 
  # as.data.frame() %>% select(Value=V1) %>% rownames_to_column('Sample') %>%
  mutate(count=h1) %>% filter(count>1) %>% #remove 23-1
  separate(Sample,c('patient','timepoint')) %>% add_clinical()
summary(lm(Value~count,h1)) #p=2e-6, adj r^2=0.35
ggplot(h1,aes(count,Value))+geom_point(aes(color=clinical))+
  scale_y_log10()+scale_x_log10()+geom_smooth(method='lm')
mutate(h1,resid=lm(Value~count,h1)$residuals) %>% 
  ggplot(aes(reorder(interaction(patient,timepoint),desc(resid)),
             resid,fill=clinical))+
  geom_col()+theme(axis.text.x=element_text(angle=90))
# same kind of patterns as TM pbmc - inconclusive overall, negative for most resid, positive on gini

# aight fuck it lazy 821 plot time
map_dfr(names(rep4$data[rep4$meta$timepoint %in% c('1')]),function(x){
  y <- na.omit(rep4$data[[x]]$expTil)
  p <- str_extract(x,'.*(?=-)')
  tp <- str_extract(x,'.$')
  data.frame('patient'=p,'timepoint'=tp,
             'high'=mean(y>=8),'low'=mean((y>1)&(y<8)),'singlet'=mean(y==1))
}) %>% add_clinical() %>% 
  pivot_longer(c(high,low,singlet),names_to='expansion') %>% 
  plotBoxplot('value')+facet_wrap(vars(expansion))
# yeah still shows a lot more high in lynch, but not as pretty
#### do like a gradient of color in a column but idk how

# 821 plot for BM expansion of BM tils
map_dfr(names(rep4$data[rep4$meta$timepoint %in% c('1')]),function(x){
  y <- rep4$data[[x]]$Clones[!is.na(rep4$data[[x]]$expTil)]
  p <- str_extract(x,'.*(?=-)')
  tp <- str_extract(x,'.$')
  data.frame('patient'=p,'timepoint'=tp,
             'high'=mean(y>=8),'low'=mean((y>1)&(y<8)),'singlet'=mean(y==1))
}) %>% add_clinical() %>% 
  pivot_longer(c(high,low,singlet),names_to='expansion') %>% 
  plotBoxplot('value')+facet_wrap(vars(expansion))
# seems like actually more clonal nr, more singlet in responder

# does the earlier regression stuff indicate lower evenness in responders?
h1 <- map_dfr(rep4$data[(rep4$meta$timepoint %in% c('1'))&
                          (!(names(rep4$data) %in% c('23-1','20-5')))],function(x){
                            x <- x[!is.na(x$expTil),]
                            if(any(duplicated(x$CDR3.aa))) warning('Duplicated CDR3aa not accounted for')
                            data.frame('div'=vegan::diversity(x$Clones,'shannon'),
                                       'nClones'=nrow(x),'nSeqs'=sum(x$Clones),
                                       'patient'=x$patient[1],'timepoint'=x$timepoint[1])
                          }) %>% add_clinical()
plotBoxplot(h1,'div')
summary(lm(div~nClones,h1))
ggplot(h1,aes(nClones,div))+geom_point(aes(color=clinical))+
  scale_y_log10()+scale_x_log10()+geom_smooth(method='lm')
mutate(h1,resid=lm(div~nClones,h1)$residuals) %>% 
  ggplot(aes(reorder(interaction(patient,timepoint),desc(resid)),
             resid,fill=clinical))+
  geom_col()+theme(axis.text.x=element_text(angle=90))
# no because there isn't even a regression difference in diversity


# predictive diversity

# tp 0/1 diversity
rep4$data[rep4$meta$timepoint=='1'] %>% 
  # map(~filter(.x,cluster_fine=='CD8 TEMRA')) %>%
  repDiversity(.method='inv.simp') %>%
  mutate(patient=str_extract(Sample,'.*(?=-)')) %>% 
  add_clinical() %>% plotBoxplot()+scale_y_log10()+
  # ggrepel::geom_text_repel(aes(label=patient),max.overlaps=25)+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1))
# lower div lynch for 1, holds for residualized version
# l>nlr>nr div for 0