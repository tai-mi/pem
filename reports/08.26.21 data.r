
# initialize --------------------------------------------------------------

knitr::opts_chunk$set(echo=F,message=F,warning=F)

library(magrittr)
library(CellaRepertorium)
library(immunarch)
library(tidyverse)

rep4 <- readRDS('data_misc/rep4.rds')
cella <- readRDS('data_misc/cella4.rds')
cl <- readRDS('data_misc/clinical.rds')
# dat <- readRDS('data_sc/meta_ldat.rds') %>% 
dat <- readRDS('data_sc/meta_tc.rds')

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
plotErrorbar2 <- function(df,value_col='Value',type=1){
  if(type==1){
    ggplot(df,aes(clinical,!!rlang::sym(value_col),
                  fill=clinical,color=clinical)) %>% prism()+
      geom_col()+
      geom_errorbar(aes(ymax=Value+serr,ymin=Value-serr),width=0.2,lwd=1)+
      theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  } else if(type==2){
    ggplot(df,aes(clinical,!!rlang::sym(value_col),
                  fill=clinical,color=clinical)) %>% prism()+
      geom_col()+
      geom_errorbar(aes(ymax=Value+serr,ymin=Value-serr),width=0.2,lwd=1)+
      theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.x=element_blank())+labs(x='')
  }
  else warning('invalid type')
}
plotErrorbar1 <- function(df,value_col='Value',plot=F,...){
  # need to pass pre-grouped DF where each row is an observation
  out <- summarize(df,serr=sd(!!rlang::sym(value_col),na.rm=T)/
                     sqrt(sum(!is.na(!!rlang::sym(value_col)))),
                   Value=mean(!!rlang::sym(value_col),na.rm=T)) %>% 
    ungroup()
  if(plot) plotErrorbar2(out,...)
  else out
}

library(ggprism)
cols_clinical <- c('Lynch-like responder'='#D25B3F','Nonlynch-like responder'='#356870','Nonresponder'='#8D7233','Healthy'='#3C367E')
cols_clinical_dark <- colorspace::darken(cols_clinical,0.6) %>% set_names(c('Lynch-like responder','Nonlynch-like responder','Nonresponder','Healthy'))
prism <- function(gg){
  gg+theme_prism()+scale_fill_manual(values=cols_clinical)+
    scale_color_manual(values=cols_clinical_dark)
}


# excel -------------------------------------------------------------------

library(openxlsx)
wb <- createWorkbook()
addSheet <- function(df,nm){
  tryCatch(addWorksheet(wb,nm),error=function(e){
    if(str_ends(e,'does not exist.')) {removeWorksheet(wb,nm);addWorksheet(wb,nm)}
    else stop(e)})
  writeData(wb,nm,df)
}


# pre/post/change in diversity --------------------------------------------

h1 <- rep4$data[!(names(rep4$data) %in% c('20-5','23-1'))] %>% 
  repDiversity('inv.simp') %>% 
  separate(Sample,c('patient','timepoint')) %>% add_clinical()
# pre
filter(h1,timepoint %in% c('1','h')) %>% 
  addSheet('pretreatment diversity')
# post
filter(h1,timepoint %in% c('3','5','h')) %>% 
  addSheet('posttreatment diversity')


# CD8 T cell prop ---------------------------------------------------------

h1 <- filter(dat,cluster_coarse!='Junk') %>% 
  group_by(patient,timepoint) %>% 
  summarize(`CD8 TEMRA`=mean(cluster_fine=='CD8 TEMRA'),
            `CD8 Effector Memory`=mean(cluster_fine=='CD8 Effector Memory')) %>% 
  add_clinical() %>% filter(interaction(patient,timepoint) %!in% c('20.5','23.1'))
# proportions of cd8 activated 
mutate(h1,value=`CD8 TEMRA`+`CD8 Effector Memory`) %>% 
  addSheet('Cell proportion - CD8 Activated')
# proportions of cd8 activated split
h1 %>% addSheet('Cell proportion - CD8 Act split')


# CD8 overall prop --------------------------------------------------------

h2 <- readRDS('data_sc/meta_ldat.rds') %>% 
  mutate(temp=paste0(patient,'-',timepoint)) %>% 
  .$temp %>% table()
h1 <- table(paste0(dat$patient,'-',dat$timepoint),dat$cluster_fine) %>% 
  as.data.frame() %>% mutate(patient=str_extract(Var1,'^.*(?=-)')) %>% 
  rename(type=Var2) %>% add_clinical() %>% filter(Var1 %!in% c('20-5','23-1')) %>% 
  rowwise() %>% mutate(value=Freq/h2[Var1])
# CD8 overall
filter(h1,type %in% c('CD8 TEMRA','CD8 Effector Memory')) %>% 
  group_by(Var1,clinical,patient) %>% summarize(value=sum(value)) %>% 
  addSheet('CD8 overall prop')
# CD8 overall split
filter(h1,type %in% c('CD8 TEMRA','CD8 Effector Memory')) %>% 
  addSheet('CD8 overall prop split')


# rarefaction -------------------------------------------------------------

clone_prop <- function(df,n,margin='seqs'){
  if(margin=='seqs') sum(filter(df,Clones>=n)$Proportion)
  else if(margin=='clones') mean(df$Clones>=n)
  else stop('invalid margin')
}
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
    summarize(serr=sd(prop)/sqrt(count[1]),prop=mean(prop)) %>% ungroup()
}
# pbmc rarefaction
clone_curve2(rep4$data[(rep4$meta$timepoint %in% c('1','3','5','h'))&
                         (!(names(rep4$data) %in% c('20-5','23-1')))]) %>% 
  addSheet('pbmc raref seqs')
# TM rarefaction
h3 <- rep4$data[(rep4$meta$timepoint %in% c(1,3,5))] %>% 
  map(~mutate(.x,Clones=expTil) %>% filter(!is.na(expTil)))
h2 <- map_dbl(h3,nrow)
h3 <- h3[h2>=20] # filter at least that many clones
map(h3,~mutate(.x,Proportion=proportions(Clones))) %>% clone_curve2(max1=100) %>% 
  addSheet('TM raref seqs')


# prop TM seqs ------------------------------------------------------------

# proportion of overall seqs
h1 <- rep4$data[(rep4$meta$timepoint %in% c(1,3,5))&
                  (!(names(rep4$data) %in% c('20-5','23-1')))]
map_dbl(h1,~sum(.x$Clones[!is.na(.x$expTil)])/sum(.x$Clones)) %>%
  set_names(names(h1)) %>% data.frame('prop'=`.`) %>%
  rownames_to_column() %>% separate(rowname,c('patient','timepoint')) %>% 
  add_clinical() %>% addSheet('TM prop overall')
# cell type proportions
group_by(dat,patient,timepoint,cluster_fine) %>% 
  summarize(prop=mean(!is.na(expTil))) %>% add_clinical() %>% 
  addSheet('TM prop cell type')


# finish ------------------------------------------------------------------

saveWorkbook(wb,'../../../../../Downloads/082621_data.xlsx',overwrite=T)
