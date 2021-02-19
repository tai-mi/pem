
# libraries ---------------------------------------------------------------

library(Seurat)
library(tidyverse)

# stuff -------------------------------------------------------------------

meta <- readRDS('sc_data/metadata.rds')
mark <- readRDS('sc_data/all_markers.rds')

mark$positivity <- ifelse(mark$avg_logFC>0,T,F)

mark$cluster %>% table() %>% view()

topMarkers <- map_dfr(levels(mark$cluster),function(l){
  mark[mark$cluster==l,][1:10,]
})

sf <- function(genes,df=mark){
  rownames(df)[str_detect(rownames(df),toupper(genes))]
}
sf1 <- function(genes,df=mark){
  for (g in genes) {
    print(rownames(df)[str_detect(rownames(df),toupper(g))])
  }
}
# ff <- function(genes,df=mark){
#   df[toupper(genes)[toupper(genes) %in% rownames(df)],]
# }
ff <- function(genes,df=mark){
  g2 <- map(genes,function(g) {
    rownames(df)[str_detect(rownames(df),toupper(g))]
  }) %>% unlist()
  df[toupper(g2)[toupper(g2) %in% rownames(df)],]
}

{
  h18 <- mark %>% filter(cluster==18)
  h27 <- mark %>% filter(cluster==27)
  g18 <- c('ITGAX','NRP1','clec4c','il3ra','lilra4','tlr7','tlr9')
  # 18 has really strong pDC markers (no ITGAX, none have tlrs), 27 has very few of these though even though very close on umap
}
{
  h16 <- mark %>% filter(cluster==16)
  g16 <- c('CD1C','CD83','CD209','THBD','CD8a','ly75','itgae','itgam','itgax')
}