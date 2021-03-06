# libraries ---------------------------------------------------------------

# library(Seurat)
library(tidyverse)

# stuff -------------------------------------------------------------------

# meta <- readRDS('data_sc/metadata.rds')
mark2 <- readRDS('data_sc/all_markers2.rds')

mark2$positivity <- ifelse(mark2$avg_logFC>0,T,F)

mark2$cluster %>% table() %>% view()

topMarkers <- map_dfr(levels(mark2$cluster),function(l){
  mark2[mark2$cluster==l,][1:10,]
})

sf <- function(genes,df=mark2){
  rownames(df)[str_detect(rownames(df),toupper(genes))]
}
sf1 <- function(genes,df=mark2){
  for (g in genes) {
    print(rownames(df)[str_detect(rownames(df),toupper(g))])
  }
}
# ff <- function(genes,df=mark2){
#   df[toupper(genes)[toupper(genes) %in% rownames(df)],]
# }
ff <- function(genes,df=mark2){
  g2 <- map(genes,function(g) {
    rownames(df)[str_detect(rownames(df),toupper(g))]
  }) %>% unlist()
  df[toupper(g2)[toupper(g2) %in% rownames(df)],]
}

{
  h18 <- mark2 %>% filter(cluster==18)
  h27 <- mark2 %>% filter(cluster==27)
  g18 <- c('ITGAX','NRP1','clec4c','il3ra','lilra4','tlr7','tlr9')
  # 18 has really strong pDC markers (no ITGAX, none have tlrs), 27 has very few of these though even though very close on umap
}
{
  h16 <- mark2 %>% filter(cluster==16)
  g16 <- c('CD1C','CD83','CD209','THBD','CD8a','ly75','itgae','itgam','itgax')
}


# Bcells ------------------------------------------------------------------

bcomp <- function(genes){
  genes <- limma::alias2Symbol(genes)
  eight <- mark2 %>% filter(cluster==8)
  rownames(eight) %<>% str_extract('^[^\\.]+')
  twosix <- mark2 %>% filter(cluster==26)
  rownames(twosix) %<>% str_extract('^[^\\.]+')
  for(g in genes){
    g <- toupper(g)
    q <- ifelse(g %in% rownames(eight),
           eight1 <- eight[g,2],
           eight1 <- NA)
    q <- ifelse(g %in% rownames(twosix),
           twosix1 <- twosix[g,2],
           twosix1 <- NA)
    cat('8:',eight1,'\n','26:',twosix1,'\n')
  }
}
# fuck it

#