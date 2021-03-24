# libraries ---------------------------------------------------------------

library(magrittr)
library(tidyverse)
# library(Seurat)


# import --------------------------------------------------------------------

mark <- map(list.files('data_sc/mark2/'), function(p){
  p2 <- paste0('data_sc/mark2/',p)
  readRDS(p2)
})
# mark %<>% set_names(paste0('c', list.files('data_sc/mark/') %>% 
                             # str_extract('\\d+')))
mark %<>% set_names(paste0('c', list.files('data_sc/mark2/') %>% 
                             str_extract('(?<=mark2_)\\d+')))


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

celldat <- data.frame(clusters=0:(length(list.files('data_sc/mark2'))-1))
# g <- c('cd14','cd19','cd3g','cd3e','cd3d','ccr7','ms4a1','il2rb','il2rg','tbx21','cd4','cd8a','cd8b','sell')
g <- cd$gene2
sf(g)
for (i in g) {
  celldat <- sf(i,return=T,xcel=T) %>% set_names(c(i,'clust')) %>% 
    inner_join(celldat,.,by=c('clusters'='clust'))
}
# removing duplicate genes
celldat %<>% .[,names(celldat)[!str_detect(names(celldat),'\\.y')]]
celldat %<>% set_names(names(celldat) %>% str_remove('\\.x'))

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
  celldat %>% column_to_rownames('clusters') %>% 
    .[,cellGenes[[ct]]] #%T>% view()
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


# panglao/rowmeans --------------------------------------------------------

dat <- readRDS('data_sc/dat_scaled2.rds')
dat %<>% as.data.frame() %>% set_colnames(paste0('c',0:33)) %>% 
  multiply_by(100)

require(magrittr)
require(tidyverse)
# require(Seurat)

db <- data.table::fread('data_markers/PanglaoDB_markers_27_Mar_2020.tsv',
                        colClasses=c(organ='factor','cell type'='factor'))

db %<>% filter(species %in% c('Hs','Mm Hs')) %>% 
  select(symbol='official gene symbol',
         type='cell type',
         ubiq='ubiquitousness index',
         description='product description',
         organ, 
         sensitivity=sensitivity_human,
         specificity=specificity_human)

map_dfc(levels(db$type),function(n){
  feats <- filter(db,type==n) %>% .$symbol
  dat[rownames(dat) %in% feats,] %>% colSums()
  # apply(dat[rownames(dat) %in% feats,], 2, median)
}) %>% set_colnames(levels(db$type)) %>% set_rownames(paste0('c',0:33)) %>% 
write.table('data_markers/sum1.csv',sep=',')

qdb <- function(gs,ordering=T){
  # returns the expression data for the genes (gs) ordered by summed expression
  # require(gt)
  gs <- toupper(gs)
  gs[!(gs %in% rownames(dat))] %>% paste('Not present:',.) %>% 
    print()
  gs <- gs[gs %in% rownames(dat)]
  if(length(gs)==1){
    if(ordering) dat[gs,order(-dat[gs,])]
    else dat[gs,]
  } else{
    if(ordering) dat[gs,order(-colSums(dat[gs,]))] #%>% gt(rownames_to_stub=T)
    else dat[gs,] #%>% gt(rownames_to_stub=T)
  }
}

#' 11,25-27,29,31,32 - B cell
#' 28? - plasma
#' 1,6,17 - NK/effector thingy?
#' 0,3,4,9,12,14?,16?,23 - t cell
#' 2 - cd8 t/NK
#' 5,8,10,15,19?,30? - DC
#' 24 - DC/pDC
#' 21 - nk, gdT, or T
#' 7,13?,20? - mt/dead/dying
#' 18 - platelets
#' 22 - dividing?
#' 33 - erythrocyte
#' 
#' CD3+ - 23x,32,16x,14x,4x,26,28,0x,10,18,8,3x,5,17,2,9x,24,15,1,12x,6,25,19,13,31,7,22,27,20,11,29,21,30,33