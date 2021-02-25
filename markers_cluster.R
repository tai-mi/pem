# libraries ---------------------------------------------------------------

library(magrittr)
library(tidyverse)
library(Seurat)


# import --------------------------------------------------------------------

mark <- map(list.files('data_sc/mark/'), function(p){
  p2 <- paste0('data_sc/mark/',p)
  readRDS(p2)
})
mark %<>% set_names(paste0('c', list.files('data_sc/mark/') %>% 
                             str_extract('\\d+')))


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
        mutate(positivity=avg_logFC>0)
      if(xcel){
        dat %>% select(avg_logFC,clust) %>% 
          full_join(tibble(clust=0:27),by=c('clust'='clust')) %>% 
          arrange(clust) %>% .$avg_logFC
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

# CD ----------------------------------------------------------------------

cd <- data.table::fread('data_markers/cd.csv')
cd$gene2 <- map_chr(1:nrow(cd),function(i){
  symbol <- limma::alias2Symbol(cd[[i,1]])
  if(identical(symbol,character(0))){
    return('NA')
  }
  else if(length(symbol)>1){
    return(symbol[1])
  }
  else if(length(symbol)==1){
    return(symbol)
  }
  else{return('NA')}
})

cd %<>% filter(gene2 != 'NA')
cd$cell2 <- map(cd$cell,function(c){
  str_split(c,',') %>% unlist() %>% str_trim()
})

# celldat -------------------------------------------------------------------

celldat <- data.frame(clusters=0:27)
# g <- c('cd14','cd19','cd3g','cd3e','cd3d','ccr7','ms4a1','il2rb','il2rg','tbx21','cd4','cd8a','cd8b','sell')
g <- cd$gene2
sf(g)
for (i in g) {
  celldat[[toupper(i)]] <- sf(i,return=T,xcel=T)
}

# need to add in data on which cell types for each gene
# probably make uhhh 

cellTypes <- unlist(cd$cell2) %>% unique() %>% sort()

cellGenes <- list()
for(c in cellTypes){
  cellGenes[[c]] <- map(1:nrow(cd),function(i){
    if(c %in% cd$cell2[[i]]) return(cd$gene2[i])
  }) %>% unlist()
}

cellFilt <- function(ct){
  if(!(ct %in% cellTypes)) break
  celldat[,cellGenes[[ct]]] %T>% view()
}
cellFilt2 <- function(ct){
  if(!(ct %in% cellTypes)) break
  celldat[,cellGenes[[ct]]] %>%
    as.matrix() %>% rowSums(na.rm=T)
}

cellScores <- list()
for(c in cellTypes){
  cellScores[[c]] <- cellFilt(c) %>% as.matrix() %>% rowSums(na.rm=T)
}
cellScores %<>% as.data.frame()
