# library1 -----------------------------------------------------------------

library(magrittr)
library(CellaRepertorium)
library(immunarch)
library(tidyverse)

rep4 <- readRDS('data_misc/rep4.rds')
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
`%!in%` <- function(lhs,rhs) !(lhs %in% rhs)

plotBoxplot <- function(df,value_col='Value'){
  df$Value <- df[[value_col]]
  ggplot(df,aes(clinical,Value,fill=clinical))+
    geom_boxplot(outlier.alpha=0)+geom_jitter(width=0.05,height=0)
}
plotErrorbar <- function(group_df,value_col='value',plot='standard'){
  group_df <- summarize(group_df,
                        serr=sd(!!rlang::sym(value_col),na.rm=T)/
                          sqrt(sum(!is.na(!!rlang::sym(value_col)))),
                        value=mean(!!rlang::sym(value_col),na.rm=T))
  if(plot=='standard'){
    ggplot(group_df,aes(clinical,value,fill=clinical,color=clinical)) %>% 
      prism()+geom_col()+
      geom_errorbar(aes(ymax=value+serr,ymin=value-serr),width=0.2,lwd=1)
  } else if(plot=='timepoint col'){
    ggplot(group_df,aes(timepoint,value,fill=clinical,color=clinical)) %>%
      prism()+geom_col(position=position_dodge(width=1))+
      geom_errorbar(aes(ymax=value+serr,ymin=value-serr),width=0.2,lwd=1,
                    position=position_dodge(width=1))
  } else if(plot=='timepoint point'){
    ggplot(group_df,aes(timepoint,value,fill=clinical,color=clinical))+
      geom_point(size=2.5,position=position_dodge(width=0.15))+
      geom_errorbar(aes(ymax=value+serr,ymin=value-serr),width=0.2,lwd=1,
                    position=position_dodge(width=0.15))+
      theme_prism()+scale_color_manual(values=cols_clinical)
  } else if(plot=='timepoint line'){
    ggplot(group_df,aes(timepoint,value,fill=clinical,color=clinical))+
      geom_point(size=2.5,position=position_dodge(width=0.15))+
      geom_line(aes(group=clinical),position=position_dodge(width=0.15),lwd=1)+
      geom_errorbar(aes(ymax=value+serr,ymin=value-serr),width=0.2,lwd=1,
                    position=position_dodge(width=0.15))+
      theme_prism()+scale_color_manual(values=cols_clinical)
  } else group_df
}

library(ggprism)
cols_clinical <- c('Lynch-like responder'='#D25B3F','Nonlynch-like responder'='#356870','Nonresponder'='#8D7233','Healthy'='#3C367E')
cols_clinical_dark <- colorspace::darken(cols_clinical,0.6) %>% set_names(c('Lynch-like responder','Nonlynch-like responder','Nonresponder','Healthy'))
prism <- function(gg){
  gg+theme_prism()+scale_fill_manual(values=cols_clinical)+
    scale_color_manual(values=cols_clinical_dark)
}


# sample overlap ----------------------------------------------------------

# repoverlap til only
repOverlap(rep4$data[rep4$meta$timepoint %in% c('0','7')],'jaccard') %>% 
  vis('heatmap2')
# repoverlap all tcr
h1 <- repOverlap(rep4$data,'overlap')
pheatmap::pheatmap(h1,
  color=colorRampPalette(c("#67001f","#d6604d","#f7f7f7",
                           "#4393c3","#053061"))(1000),
  breaks=c(seq(min(h1,na.rm=T),0.09, length.out=1001),max(h1,na.rm=T)),
  annotation_col=data.frame('clinical'=factor(rep4$meta$clinical),
                            row.names=names(rep4$data)))


# cell types of clonotypes ------------------------------------------------

#### modify to add partial counts when stuff tied for primary
nfilter <- 3 #only clones with at least this many IDed cells
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
  ggplot(aes(primary,value,fill=secondary))+geom_col(position='stack')
# don't have scaling to correct for the number of cells of each type yet
# but treg lot of overlap cd4act cd8, cd8exh-act strong, cd4 act-exh?
# lot more treg and cd4 act in nlr, low treg in l, 


# clonotype size rarefaction ----------------------------------------------

clone_curve2 <- function(dfs,margin='seqs',max1=20){
  clone_prop <- function(df,n,margin='seqs'){
    if(margin=='seqs') sum(filter(df,Clones>=n)$Proportion)
    else if(margin=='clones') mean(df$Clones>=n)
    else stop('invalid margin')
  }
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


# proportion seqs TM temporal change ----------------------------------------

h1 <- rep4$data[(rep4$meta$timepoint %in% c(1,3,5))&
                  (!(names(rep4$data) %in% c('20-5','23-1')))]
h1 <- map_dbl(h1,~sum(.x$Clones[!is.na(.x$expTil)])/sum(.x$Clones)) %>%
  set_names(names(h1)) %>% data.frame('prop'=`.`) %>%
  rownames_to_column() %>% separate(rowname,c('patient','timepoint')) %>% 
  add_clinical()
h1 <- pivot_wider(h1,id_cols=patient,names_from=timepoint,values_from=prop) %>% 
  mutate(`3-1`=(`3`-`1`)/`1`,`5-1`=(`5`-`1`)/`1`,`5-3`=(`5`-`3`)/`3`) %>% 
  add_clinical()
plotBoxplot(h1,'5-1')
filter(h1,!is.na(`5-1`)) %>% 
  ggplot(aes(reorder(patient,desc(`5-1`)),`5-1`,fill=clinical))+geom_col()


# proportion TM seqs by cell type -----------------------------------------

group_by(dat,patient,timepoint,cluster_fine) %>% 
  summarize(prop=mean(!is.na(expTil))) %>% add_clinical() %>% ungroup() %>% 
  group_by(clinical,cluster_fine) %>% 
  summarize(value=mean(prop),serr=sd(prop)/sqrt(n())) %>% 
  filter(clinical!='Healthy') %>% 
  ggplot(aes(clinical,value,fill=clinical,color=clinical)) %>% 
  prism()+geom_col()+
  geom_errorbar(aes(ymax=value+serr,ymin=value-serr),width=0.2,lwd=1)+
  facet_wrap(~cluster_fine,scales='free')+
  theme(axis.text.x=element_blank())
# coarse
group_by(dat,patient,timepoint,cluster_coarse) %>% 
  summarize(prop=mean(!is.na(expTil))) %>% add_clinical() %>% ungroup() %>% 
  group_by(clinical,cluster_coarse) %>% 
  summarize(value=mean(prop),serr=sd(prop)/sqrt(n())) %>% 
  filter(clinical!='Healthy') %>% 
  ggplot(aes(clinical,value,fill=clinical,color=clinical)) %>% 
  prism()+geom_col()+
  geom_errorbar(aes(ymax=value+serr,ymin=value-serr),width=0.2,lwd=1)+
  facet_wrap(~cluster_coarse,scales='free')+
  theme(axis.text.x=element_blank())

# coarse
group_by(dat,patient,timepoint,cluster_coarse) %>% 
  summarize(prop=mean(!is.na(expTil))) %>% add_clinical() %>% ungroup() %>% 
  group_by(clinical,timepoint,cluster_coarse) %>% 
  summarize(value=mean(prop),serr=sd(prop)/sqrt(n())) %>% 
  filter(clinical!='Healthy') %>% 
  ggplot(aes(timepoint,value,fill=clinical,color=clinical))+
  geom_point(size=2.5,position=position_dodge(width=0.15))+
  geom_line(aes(group=clinical),position=position_dodge(width=0.15),lwd=1)+
  geom_errorbar(aes(ymax=value+serr,ymin=value-serr),width=0.2,lwd=1,
                position=position_dodge(width=0.15))+
  theme_prism()+scale_color_manual(values=cols_clinical)+
  facet_wrap(~cluster_coarse,scales='free')+
  theme(axis.text.x=element_blank())
# fine
group_by(dat,patient,timepoint,cluster_fine) %>% 
  summarize(prop=mean(!is.na(expTil))) %>% add_clinical() %>% ungroup() %>% 
  group_by(clinical,timepoint,cluster_fine) %>% 
  summarize(value=mean(prop),serr=sd(prop)/sqrt(n())) %>% 
  filter(clinical!='Healthy') %>% 
  ggplot(aes(timepoint,value,fill=clinical,color=clinical))+
  geom_point(size=2.5,position=position_dodge(width=0.15))+
  geom_line(aes(group=clinical),position=position_dodge(width=0.15),lwd=1)+
  geom_errorbar(aes(ymax=value+serr,ymin=value-serr),width=0.2,lwd=1,
                position=position_dodge(width=0.15))+
  theme_prism()+scale_color_manual(values=cols_clinical)+
  facet_wrap(~cluster_fine,scales='free')+
  theme(axis.text.x=element_blank())
# pre vs post
mutate(dat,timepoint=if_else(timepoint=='1','pre','post') %>%
         factor(levels=c('pre','post'))) %>% 
  group_by(patient,timepoint,cluster_fine) %>% 
  summarize(prop=mean(!is.na(expTil))) %>% add_clinical() %>% ungroup() %>% 
  group_by(clinical,timepoint,cluster_fine) %>% 
  summarize(value=mean(prop),serr=sd(prop)/sqrt(n())) %>% 
  filter(clinical!='Healthy') %>% 
  ggplot(aes(timepoint,value,fill=clinical,color=clinical))+
  geom_point(size=2.5,position=position_dodge(width=0.15))+
  geom_line(aes(group=clinical),position=position_dodge(width=0.15),lwd=1)+
  geom_errorbar(aes(ymax=value+serr,ymin=value-serr),width=0.2,lwd=1,
                position=position_dodge(width=0.15))+
  theme_prism()+scale_color_manual(values=cols_clinical)+
  facet_wrap(~cluster_fine,scales='free')+
  theme(axis.text.x=element_blank())
