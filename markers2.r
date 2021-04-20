# libraries ---------------------------------------------------------------

library(magrittr)
library(tidyverse)

# stuff -------------------------------------------------------------------

allMark <- readRDS('data_sc/tc_allmark.rds')
topMark <- map(group_split(allMark,cluster),~arrange(.x,desc(avg_log2FC)) %>% .[1,]) %>% 
  set_names(group_by(allMark,cluster) %>% group_keys() %>% unlist()) %>% 
  map_dfr(~.x)

source('panglao.r')
panglao2 <- function(clust){
  allMark %>% filter(cluster==clust) %>% arrange(desc(avg_log2FC)) %>% 
    .$gene %>% .[1:10] %>% 
    panglaoMain()
}


# stuff2 ------------------------------------------------------------------

mark <- map(list.files('data_sc/marktc/'), function(p){
  p2 <- paste0('data_sc/marktc/',p)
  readRDS(p2)
})
# mark %<>% set_names(paste0('c', list.files('data_sc/mark/') %>% 
# str_extract('\\d+')))
mark %<>% set_names(paste0('c', list.files('data_sc/marktc/') %>% 
                             str_extract('(?<=mark_)\\d+')))


# marker functions --------------------------------------------------------

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

