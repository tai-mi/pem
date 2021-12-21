# libraries
library(magrittr)
library(CellaRepertorium)
library(immunarch)
library(tidyverse)

# base functions
cl <- readRDS('data_misc/clinical.rds')
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
`%!in%` <- function(lhs,rhs) !(lhs %in% rhs)

# load functions
load_tcr <- function(){
  rep4 <- readRDS('data_misc/rep4.rds')
  rep4$data <- rep4$data[rep4$meta$patient!='24']
  rep4$meta <- filter(rep4$meta,patient!='24')
  rep4$meta$clinical[rep4$meta$timepoint=='h'] <- 'Healthy'
  assign('rep4',rep4,envir=.GlobalEnv)
  
  cella <- readRDS('data_misc/cella4.rds')
  cella@contig_tbl %<>% filter(patient!='24')
  cella@cell_tbl %<>% filter(patient!='24')
  cella@contig_tbl %<>% add_clinical()
  cella@cell_tbl %<>% add_clinical()
  assign('cella',cella,envir=.GlobalEnv)
}
load_main <- function(type){
  if(type=='ldat') dat <- readRDS('data_sc/meta_ldat.rds')
  if (type=='tc') dat <- readRDS('data_sc/meta_tc.rds') %>% filter(patient!='24')
  dat %<>% filter(patient!='24')
}

# plotting functions
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
                  fill=clinical,color=clinical)) %>% prism()+geom_col()+
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
  if(plot) plotErrorbar2(out,value_col='Value',...)
  else out
}
cols_clinical <- c('Lynch-like responder'='#D25B3F','Nonlynch-like responder'='#356870','Nonresponder'='#8D7233','Healthy'='#3C367E')
cols_clinical_dark <- colorspace::darken(cols_clinical,0.6) %>% set_names(c('Lynch-like responder','Nonlynch-like responder','Nonresponder','Healthy'))
prism <- function(gg){
  require(ggprism)
  gg+theme_prism()+scale_fill_manual(values=cols_clinical)+
    scale_color_manual(values=cols_clinical_dark)
}