# libraries ---------------------------------------------------------------

library(magrittr)
library(tidyverse)

# stuff -------------------------------------------------------------------

allMark <- readRDS('data_sc/tc_allmark.rds')
topMark <- map(group_split(allMark,cluster),~arrange(.x,desc(avg_log2FC)) %>% .[1,]) %>% 
  set_names(group_by(allMark,cluster) %>% group_keys() %>% unlist()) %>% 
  map_dfr(~.x)

source('../rig-i/panglao.r')
panglao2 <- function(clust){
  allMark %>% filter(cluster==clust) %>% arrange(desc(avg_log2FC)) %>% 
    .$gene %>% .[1:10] %>% 
    panglaoMain()
}


# mark ------------------------------------------------------------------

mark <- map(list.files('data_sc/mark_tc/'), function(p){
  p2 <- paste0('data_sc/mark_tc/',p)
  readRDS(p2)
})
# mark %<>% set_names(paste0('c', list.files('data_sc/mark/') %>% 
# str_extract('\\d+')))
mark %<>% set_names(paste0('c', list.files('data_sc/mark_tc/') %>% 
                             str_extract('(?<=mark_)\\d+')))


# marker functions --------------------------------------------------------

sf2 <- function(gene){
  gene <- toupper(gene)
  colnames(scaled)[str_detect(colnames(scaled),gene)]
}

sf <- function(genes,id=NULL,return=F,xcel=F){
  if(!is.null(id)){
    df <- mark[[id]]
    if(return){
      mark[[id]][rownames(mark[[id]]) %in% toupper(genes),]
    }
    else{
      for(g in genes){
        print(rownames(df)[str_detect(rownames(df),toupper(g))])
      }
    } 
  }
  else{
    if(return){
      dat <- map_dfr(names(mark),function(n){
        m <- mark[[n]]
        m[rownames(m) %in% toupper(genes),] %>% as.data.frame() %>% 
          mutate(clust=str_extract(n,'[0-9]+') %>% as.integer())
      }) %>% 
        mutate(positivity=avg_log2FC>0)
      if(xcel){
        dat %>% select(avg_log2FC,clust) %>% 
          full_join(tibble(clust=0:(length(list.files('data_sc/mark2'))-1)),by=c('clust'='clust')) #%>% 
        # arrange(clust) %>% .$avg_logFC
        # paste(.,'',collapse='')
      }
      else{dat}
    }
    else{
      for(n in names(mark)){
        m <- mark[[n]]
        print(paste('####',n,'####'))
        for (g in genes) {
          print(rownames(m)[str_detect(rownames(m),toupper(g))])
        }
      }
    }
  }
}

dev <- function(gene){
  gene %<>% toupper()
  x <- map(mark,~.x[gene,'avg_log2FC']) %>% 
    unlist()
  x <- ifelse(is.na(x),0,x)
  x <- as_tibble(x,rownames='cluster')
  y <- mutate(x,c2=str_extract(cluster,'\\d+') %>% as.numeric()) %>% 
    # mutate(cluster=str_extract(value,'\\d+') %>% as.numeric) %>% 
    ggplot()+geom_col(aes(reorder(cluster,c2),value,fill=cluster))+
    theme(axis.text.x=element_text(angle=90,),legend.position='none')+
    ylim(c(min(c(-0.5,x$value)),max(c(0.5,x$value))))+
    geom_hline(yintercept=0)
  y
}

# expression --------------------------------------------------------------

scaled <- readRDS('data_sc/tc_scaled.rds') %>% 
  t() %>% as.data.frame()

dev2 <- function(genes,center=0,drop1=NULL){
  genes <- toupper(genes)
  if(length(genes)==1){
    x <- tibble('cluster'=as.numeric(rownames(scaled)),
           expression=scaled[[genes]])
  } else{
    x <- tibble('cluster'=as.numeric(rownames(scaled)),
           expression=scaled[,genes] %>% rowSums())
  }
  x$expression[x$cluster %in% as.character(drop1)] <- NA
  x %>% mutate(expression=expression-sum(range(expression,na.rm=T))/2-center) %>% 
    arrange(cluster) %>% 
    ggplot()+geom_col(aes(cluster,expression,fill=as.factor(cluster)))+
    theme(axis.text.x=element_text(angle=90,),legend.position='none')+
    scale_x_continuous(breaks=0:40)+
    geom_hline(yintercept=0)+
    labs(title=genes)
}

# naive: SELL, CCR7; TCF7, TIM3?; CD44-, CD69- (earlier than 44)


# eric --------------------------------------------------------------------

id1 <- readxl::read_xlsx('data_markers/tc_eric.xlsx','Sheet1') %>% 
  mutate(genes=toupper(genes) %>% replace_na(''),
         type=str_replace_all(type,'ï','i'))

addGene <- function(genes,clusts){
  genes %<>% toupper()
  if(length(genes)==1){
    for(i in 1:nrow(id1)){
      if(id1$cluster[i] %in% clusts){
        id1$genes[i] <- paste0(id1$genes[i],',',genes)
      }
    }
  } else if(length(clusts)==1){
    id1$genes[id1$cluster==clusts] <- id1[id1$cluster==clusts,'genes'] %>% 
      unlist() %>% append(genes)
  }
  id1
}

# using housekeeping genes c('gapdh','sdha','hprt1','hbs1l','ahsp','b2m'),
# 26 is really high, then 27 and 18, then kinda 14/15/28
# 6/7/9 really low, kinda 17

# https://github.com/EDePasquale/DoubletDecon