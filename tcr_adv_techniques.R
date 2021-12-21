# libraries ---------------------------------------------------------------

library(magrittr)
library(immunarch)
library(CellaRepertorium)
library(tidyverse)


# load --------------------------------------------------------------------

rep <- readRDS('data_misc/rep2_healthy.rds')
dat <- readRDS('data_misc/dat_healthy.rds')

rep$meta %<>% mutate(coarse=str_replace(timepoint,'1|3|5','pbmc'))
# rep$data <- map(rep$data,function(df){ #fix Vname issues
#   mutate(df,
#          V.name=case_when(str_detect(V.name,'/DV')~str_remove(V.name,'/'),
#                   str_detect(V.name,'/[0-9]')~str_remove(V.name,'/.*'),
#                   T~V.name))
# }) # still a bunch only in bulk; bulk has no TRAJ
# saveRDS(rep,'data_misc/rep2_healthy.rds')


# immunarch ---------------------------------------------------------------

# V gene overlap
trv <- geneUsage(rep$data,.gene='hs.trbv',.quant='count',.ambig='exc',
                  .norm=T)
# trv <- geneUsage(rep$data,.gene='hs.trav',.quant='count',.ambig='exc',
#                   .norm=T)

geneUsageAnalysis(trv,.method='tsne') %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  ggplot(aes(DimI,DimII,label=rowname,color=rep$meta$patient))+geom_point()+
  ggrepel::geom_text_repel()
geneUsageAnalysis(trv,.method='pca+hclust') %>% vis()
geneUsageAnalysis(trv,.method='js') %>% vis(.plot='heatmap2')
geneUsageAnalysis(trv,.method='js')[c(20:32,35:48,50:51,53:82),c(20:32,35:48,50:51,53:82)] %>% pheatmap::pheatmap()
geneUsageAnalysis(trv,.method='mds+dbscan')$data %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  ggplot(aes(DimI,DimII,label=rowname,color=rep$meta$coarse))+geom_point()+
  ggrepel::geom_text_repel()

# rep overlap
jc <- repOverlap(rep$data,'jaccard')
mor <- repOverlap(rep$data,'morisita')

vis(jc,.plot='heatmap2')
vis(mor,.plot='heatmap2')

repOverlapAnalysis(jc,'tsne') %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  ggplot(aes(DimI,DimII,label=rowname,color=rep$meta$coarse))+geom_point()+
  ggrepel::geom_text_repel()
repOverlapAnalysis(mor,'tsne') %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  ggplot(aes(DimI,DimII,label=rowname,color=rep$meta$lynch))+geom_point()+
  ggrepel::geom_text_repel()
repOverlapAnalysis(mor,'tsne+kmeans',.k=6) %>% vis()
#### try aggregating patients and then running this analysis on them as wholes

vdj <- data.table::fread('data_rep/vdjdb_single.tsv')
#### do stuff

# kmers
k5 <- getKmers(rep$data,5)
# k8 <- getKmers(rep$data,8)

vis(k5,.head=30,.position='fill')
kmer_profile(k5[[1]],.method='freq') %>% vis()
kmer_profile(k5[[1]],.method='freq') %>% vis(.plot='seq')


# circos ------------------------------------------------------------------

#circos repoverlap
library(circlize)
circosOverlap <- function(p){
  circos.clear()
  repOv <- rep$data[rep$meta$patient %in% p] %>% 
    repOverlap(.method='jaccard',.col='aa')
  colorFunct <- colorRamp2(
    c(max(repOv,na.rm=T),median(repOv,na.rm=T),min(repOv,na.rm=T)) %>% rev(),
    RColorBrewer::brewer.pal(3,'YlOrRd'),transparency=0.4)
  chordDiagram(repOv,col=colorFunct,directional=1,link.arr.type='big.arrow',
               direction.type=c('diffHeight','arrows'))
}
barOverlap <- function(p){ #barplot of jaccard relative to til
    repOverlap(rep$data[rep$meta$patient %in% p],
               .method='jaccard',.col='aa')[,paste0(p,'-til')] %>% 
    na.omit() %>% as.data.frame() %>% rownames_to_column() %>%
    set_colnames(c('Timepoint','Jaccard Similarity')) %>% 
    ggplot(aes(Timepoint,`Jaccard Similarity`,fill=Timepoint))+
    geom_col()
}
# barplot of jaccard relative to til for all patients
repOverlap(rep$data[(rep$meta$timepoint!='healthy')&(rep$meta$patient %in% filter(rep$meta,timepoint=='til')$patient)],'jaccard','aa') %>% 
  reshape2::melt() %>% 
  filter(str_ends(Var1,'til'),str_ends(Var2,'til',negate=T),
         str_extract(Var1,'^[0-9]+')==str_extract(Var2,'^[0-9]+')) %>% 
  mutate(timepoint=str_extract(Var2,'.$'),
         Var1=str_extract(Var1,'^[0-9]+')) %>% 
  ggplot(aes(Var1,value,fill=timepoint))+geom_col(position='dodge')+
  theme(axis.text.x=element_text(angle=-90))+
  labs(y='Jaccard Similarity',x='Patient',
       title='TIL-PBMC repertoire similarity')

#circos spectratype
circosSpectra <- function(mask1,ab){
  # mask for which of repdata and whether a or b
  require(circlize)
  vj <- map_dfr(rep$data[mask1],~.x) %>% select(Clones,J.name,V.name) %>% 
    filter(str_starts(V.name,paste0('TR',toupper(ab),'V')),
           str_starts(J.name,paste0('TR',toupper(ab),'J'))) %>% 
    na.omit() %>% group_by(J.name,V.name) %>% 
    summarise(total=sum(Clones))
  temp1 <- vj %>% group_by(V.name) %>% summarise(total=sum(total))
  vj$V.name[vj$V.name %in% temp1$V.name[temp1$total<(sum(temp1$total)/100)]] <- paste0('TR',toupper(ab),'V-other')
  temp1 <- vj %>% group_by(J.name) %>% summarise(total=sum(total))
  vj$J.name[vj$J.name %in% temp1$J.name[temp1$total<(sum(temp1$total)/100)]] <- paste0('TR',toupper(ab),'J-other')
  circos.clear()
  circos.par(gap.degree=0,
             gap.after=c(rep(0.2,length(unique(vj$J.name))-1),5,
                         rep(0.2,length(unique(vj$V.name))-1),5))
  colorFunct <- colorRamp2(
    c(max(vj$total),median(vj$total),min(vj$total)) %>% rev(),
    RColorBrewer::brewer.pal(3,'YlOrRd'),transparency=0.4)
  chordDiagram(vj,annotationTrack='grid',
               preAllocateTracks=1,
               col=colorFunct,link.sort=T,)
  circos.trackPlotRegion(track.index = 1,panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter,CELL_META$ylim[1],CELL_META$sector.index,
                facing = "clockwise", niceFacing=T,
                cex=0.7,track.index=1,adj=c(0,0.5))
    # circos.axis(h = "top", labels.cex = 0.5, major.tick.length = 0.1,
    #             sector.index=CELL_META$sector.index, track.index = 2)
  }, bg.border = NA)
}
circosSpectra(rep$meta$lynch,'b')


# cella -------------------------------------------------------------------

h1 <- cdhit_ccdb(dat,'cdr3',type='AA',cluster_pk='cdhit',
                 identity=0.80,kmerSize=5)
# h1 <- fine_clustering(dat,'cdr3',type='AA',
#                       big_memory_brute=T,method='levenshtein')

