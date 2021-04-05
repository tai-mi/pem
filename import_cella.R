# library -----------------------------------------------------------------

library(magrittr)
library(tidyverse)
library(CellaRepertorium)

# WORKFLOW - load, QC, filter, rerun QC, ...


# Load ---------------------------------------------------------------

fList <- list.files('data_tcr')
tempData <- data.frame()

for (file in fList) {
  tempData <- data.table::fread(paste0('data_tcr/',file)) %>%
    select(-one_of('V1')) %>% 
    mutate(patient=str_extract(file,'(?<=^PEM)\\d+'), 
           timepoint=str_extract(file,'(?<=C)\\d')) %>% 
    rbind(tempData,.)
}

dat <- ContigCellDB_10XVDJ(tempData,
                           contig_pk=c('barcode','contig_id','patient','timepoint'),
                           cell_pk=c('barcode','patient','timepoint'))
# ok this is lazy and has a bunch of cell characteristics in the contig tbl and duplicated id info, but it would be a mess to differentiate barcodes and transfer stuff to cell_tbl so I'm not going to

dat %<>% mutate_cdb(productive=case_when(
  toupper(as.character(productive))=='TRUE'~T,
  toupper(as.character(productive))=='FALSE'~F,
  toupper(as.character(productive))=='NONE'~NA,
))

cl <- readRDS('data_misc/clinical.rds')
clResponse <- set_names(cl$response,cl$sample)
clLynch <- set_names(cl$lynch,cl$sample)

# QC ----------------------------------------------------------------------

nrow(dat$contig_tbl)

print(paste(
  ifelse(all(dat$contig_tbl$high_confidence),
         'All high confidence,',
         paste0(round(100*mean(dat@contig_tbl$high_confidence),2),'% of cells are high confidence,')),
  ifelse(all(dat$contig_tbl$is_cell),
         'all from a cell,',paste0(round(100*mean(dat@contig_tbl$is_cell),2),'% of contigs are from cells,')),
  ifelse(all(dat$contig_tbl$full_length),
         'all full length,',
         paste0(round(100*mean(dat@contig_tbl$full_length),2),
                '% of contigs are full length,')),'and',
  ifelse(all(dat$contig_tbl$productive),'all productive',
         paste0(round(100*sum(dat@contig_tbl$productive==T,na.rm=T)
                      /nrow(dat@contig_tbl),2),'% are productive'))
))

# ggplot(dat@contig_tbl)+geom_density(aes(length))
ggplot(dat@contig_tbl)+geom_density(aes(map_dbl(cdr3,nchar)))+
  geom_vline(xintercept=8)
# slow w/ jitter (but is useful for those w/ low counts)
# ggplot(dat@contig_tbl,aes(x='reads',y=reads))+
#   geom_violin()+scale_y_log10()+#geom_jitter(alpha=0.05)+
#   facet_wrap(~patient+timepoint)
# ggplot(dat@contig_tbl,aes(x='reads',y=umis))+
#   geom_violin()+scale_y_log10()+#geom_jitter(alpha=0.05)+
#   facet_wrap(~patient+timepoint)
ggplot(dat@contig_tbl)+geom_bar(aes(x=interaction(patient,timepoint)))+
  theme(axis.text.x.bottom=element_text(angle=90))

guess_celltype(dat@contig_tbl$chain) %>% table()

# Tab_umi = crosstab_by_celltype(cdb)
# ggplot(Tab_umi, aes(color = factor(is_cell), x = T_ab, group = interaction(is_cell, sample, pop))) + stat_ecdf() + coord_cartesian(xlim = c(0, 10)) + ylab('Fraction of barcodes') + theme_minimal() + scale_color_discrete('10X called cell?')

#' Notes: 
#' very low reads in 23-15 and low in 20-55 after filter and 21-15 looks a bit light on umis -- yeah they only have 2, 65, and 270 original rows respectively
#' when it has None in the productive col, it has None for cdr3


# filter -------------------------------------------------------------------

# filter all cells, high_conf, not multichain (too many chains for one cell?), productive, full_length (?), and cdr3 > 7 aa's (also excludes Nones)
dat %<>% filter_cdb(is_cell,high_confidence,chain %in% c('TRA','TRB')
                               ,productive,full_length,nchar(cdr3)>6)

tab_umi <- crosstab_by_celltype(dat) %>% 
  full_join(dat@contig_tbl,by=c('barcode','patient','timepoint')) %>% 
  filter(!is.na(is_cell))
# adds in number of t_ab sequences for each cell barcode

saveRDS(dat,'data_misc/dat.rds')

# diversity -----------------------------------------------------------------

require(vegan)
# create diversity dataframe
divData <- data.frame(group_by(dat@contig_tbl,patient,timepoint) %>% 
                        group_keys(),
    'contigs'=group_split(dat@contig_tbl,patient,timepoint) %>% map_dbl(nrow)) %>% mutate(lynch=clLynch[patient],response=clResponse[patient])
# add shannon, invsimp and hill indices
divData$shannon <- dat@contig_tbl %>% group_split(patient,timepoint) %>% 
  map_dbl(~.x$cdr3 %>% table() %>% unlist() %>% diversity)
divData$invSimp <- dat@contig_tbl %>% group_split(patient,timepoint) %>% 
  map_dbl(~.x$cdr3 %>% table() %>% unlist() %>% diversity(index='invsimpson'))

# response %in% c('CR','PR')
# response %in% c('SD','PD')
divData %>% filter(response %in% c('SD','PD')) %>% 
  ggplot()+geom_col(aes(x=patient,y=invSimp,fill=timepoint,group=interaction(patient,timepoint)),position='dodge')+theme(axis.text.x=element_text(angle=90))+scale_y_sqrt()+labs(title='Non-Responders')

hillNums <- c(2,4,8,16)
divData$hill <- dat@contig_tbl %>% group_split(patient,timepoint) %>% 
  map(~.x$cdr3 %>% table() %>% unlist() %>% 
        renyi(scales=hillNums,hill=T))
# plot hill numbers
hillPlot <- divData$hill %>% as.data.frame() %>% t() %>%
  set_colnames(paste0(hillNums)) %>% 
  as.data.frame() %>% 
  mutate(p=divData$patient,t=divData$timepoint) %>% 
  mutate(r=clResponse[as.character(p)],l=clLynch[as.character(p)]) %>% 
  pivot_longer(cols=paste0(hillNums)) %>% 
  mutate(name=factor(name,levels=paste0(hillNums)))
hillPlot %>% filter(r %in% c('SD','PD')) %>% 
  ggplot()+geom_line(aes(name,value,col=p,group=interaction(p,t),linetype=t))+facet_wrap(~p)+scale_linetype_manual(values=c('solid','longdash','dotdash'))+labs(x='hill numbers',title='Non-Responders')+scale_y_log10()

# resample from larger numbers of contigs for the smallest number of contigs timepoint
# permutation test pairs?


# track clonotypes --------------------------------------------------------

require(streamgraph)

# counts of each cdr3 by patient x timepoint
h2 <- dat@contig_tbl %>% dplyr::count(patient,timepoint,cdr3,sort=T)
topNum <- 7 #top clones from each timepoint to look at

h3 <- h2 %>% group_split(patient) %>% map(function(x){
  tpx <- unique(x$timepoint)
  # find top `topNum` clonotypes from each timepoint for each patient
  cdr <- map(tpx,function(tp){
    y <- filter(x,timepoint==tp)
    if(nrow(y) >= topNum){
      cdr <- slice_max(y,order_by=n,n=topNum,with_ties=F) %>% .$cdr3
    }
    else cdr <- y$cdr3
  }) %>% unlist() %>% unique()
  # list mapping timepoint to number of contigs
  tpnMap <- map_dbl(tpx,
              ~filter(x,timepoint==.x) %>% .$n %>% sum()) %>% 
    set_names(as.character(tpx))
  # find all rows with top clonotypes and add proportion column
  map_dfr(cdr, ~filter(x,cdr3==.x)) %>%
    mutate(prop=100*n/(tpnMap[as.character(timepoint)]))
}) %>% 
  set_names(group_by(h2,patient) %>% group_keys() %>% unlist)

trackClones <- function(df){
  df %>% mutate(year=as.Date.character(timepoint,format='%d')) %>% 
    streamgraph('cdr3','prop','year') %>% 
    sg_axis_x(2,'day','%d') %>% 
    sg_legend(show=T)
}

for (i in clResponse[clResponse %in% c('SD','PD')] %>% names() %>% .[1:9]) {
  trackClones(h3[[i]]) %>% htmlwidgets::saveWidget('temp.html')
  webshot::webshot('temp.html',paste0('figures/nonresp/','c',i,'.png'))
}


# til connection ----------------------------------------------------------

rep1 <- readRDS('data_misc/rep1.rds')
dat <- readRDS('data_misc/dat.rds')

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

trackCl2 <- function(p,...){
  require(immunarch)
  if(!exists('rep1')) break
  p <- as.character(p)
  temp <- filter(dat@contig_tbl,patient==p) %>%
    group_split(timepoint) %>%
    map(~toImmunarch(.x))
  temp %<>% set_names(filter(dat@contig_tbl,patient==p) %>%
                        group_by(timepoint) %>% group_keys() %>% unlist())
  if(p %in% names(rep1$data)) temp$til <- rep1$data[[p]]
  trackClonotypes(temp,...)
}

# quantify proportion/prop of top clones from til that present in pbmc
# correlations between them?