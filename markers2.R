# object ------------------------------------------------------------------

library(tidyverse)

h1 <- list(TCELL = 'CD3D',

#HELPER T CELLS
TH = 'CD4',
TH1 = 'TNF,IL2,IFNG,CXCR3',
TH2 = 'IL4,IL5,IL9,IL10,CCR4',
TH17 = 'IL21,IL22,IL25,IL17A,CCR6',
TH22 = 'IL22,CCR10',
TH9 = 'IL9',
TFH = 'IL21',

#CYTOTOXIC T CELLS
TC = 'TRAC,TRBV1,TRBV2,TRBC2,CD8A,CD8B,CD28',
TC1 = 'TNF,IFNG,IL2,CXCR3,TBX21',
TC2 = 'IL4,IL5,CCR4,GATA3',
TC9 = 'IL9,IL10,IRF4',
TC17 = 'CCR6,KLRB1,IL17A,IRF4,RORC',

#REGULATOR T CELLS (cd8 or cd4)
TREG = 'FOXP3,IL2RA,STAT5A,CTLA4,IL10,TGFB1',

TNAIVE = 'ITGAE',
TEM = 'GZMA',
TRM = 'CD69,ITGAE,CTLA4',
TSCM = 'SELL',
TEFF = 'GZMA',

TNKT = 'NCAM1,FCGR3A',
TGD = 'TRGC1,TRGC2,TRDC',

#MACRO
MACROPHAGE = 'CD14,FCGR3A,FCGR1A,CD68,TFRC,CCR5',
# MORE AT "https://www.bio-rad-antibodies.com/macrophage-m1-m2-tam-tcr-FCGR3A9-cd-markers-antibodies.html"
NEUTROPHIL = 'FUT4,FCGR3A,FCGR3A,CEACAM8,CEACAM8,ITGAM,CD33,MPO',
BASOPHIL = 'CCR3,ITGAX,CD22,CD22,CD69,PTGDR2,ENPP3,FCER1A,IL3RA,ITGA2',
EOSINOPHIL = 'CCR3,ITGAM,FUT4,PTPRC,CEACAM8,ADGRE1,IL5RA,IL5RA,ITGA4,ITGA4,SIGLEC8',
# MONO
MONOCLASS = 'CD14,CCR2,CCR5,SELL',
MONONCLASS = 'CD14,FCGR3A,CX3CR1',
MONOINT = 'CD14,FCGR3A,CD68,ITGAX',
NK = 'NCAM1',
GRANULOCYTE = 'CEACAM8',
MAST = 'PTPRC,KIT,ENPP3,FCGR2A,FCGR2B,FCER1A,MITF',

#DC
DCC = 'CD1C,CD8B,CD8A,ITGAM,ITGAX,ITGAE,LY75',
DCP ='CLEC4C,NRP1,IL3RA,LILRA4,TLR7,TLR9',

#B CELL
B = 'CD19,CD24,CD40,CD72',
B1 = 'MS4A1,CD27,SPN',
BNAIVE = 'CD22',
BPLASMA = 'CD27,CD38,PRDM1,XBP1,IRF4,SDC1,BCL2,BCL6,SDC1',
BMEM = 'CD27,BCL2,NT5E',
BREG = 'IL10,CD1D,CD24,CD5,CD44,CD38,TFRC,HAVCR1')

for(i in h1){
  stringr::str_split(i,',')
}
h2 <- purrr::map(h1,function(i){
  strsplit(i,',')
})
saveRDS(h2,'searchF2.rds')



# object2 -----------------------------------------------------------------

# higher specificity
h2 <- list(TCELL = 'CD3D,CD3G,CD3E',
           
           #HELPER T CELLS
           TH = 'CD4',
           # TH1 = 'TNF,IL2,IFNG,CXCR3',
           # TH2 = 'IL4,IL5,IL9,IL10,CCR4',
           # TH17 = 'IL21,IL22,IL25,IL17A,CCR6',
           # TH22 = 'IL22,CCR10',
           # TH9 = 'IL9',
           # TFH = 'IL21',
           # 
           # #CYTOTOXIC T CELLS
           TC = 'TRAC,TRBV1,TRBV2,TRBC2,CD8A,CD8B,CD28',
           # TC1 = 'TNF,IFNG,IL2,CXCR3,TBX21',
           # TC2 = 'IL4,IL5,CCR4,GATA3',
           # TC9 = 'IL9,IL10,IRF4',
           # TC17 = 'CCR6,KLRB1,IL17A,IRF4,RORC',
           # 
           # #REGULATOR T CELLS (cd8 or cd4)
           # TREG = 'FOXP3,IL2RA,STAT5A,CTLA4,IL10,TGFB1',
           # 
           # TNAIVE = 'ITGAE',
           # TEM = 'GZMA',
           # TRM = 'CD69,ITGAE,CTLA4',
           # TSCM = 'SELL',
           # TEFF = 'GZMA',
           # 
           # TNKT = 'NCAM1,FCGR3A',
           # TGD = 'TRGC1,TRGC2,TRDC',
           # 
           # #MACRO
           MACRO = 'CD163,CD68,ITGAM',
           # # MORE AT "https://www.bio-rad-antibodies.com/macrophage-m1-m2-tam-tcr-FCGR3A9-cd-markers-antibodies.html"
           # # MONO
           MONO = 'CD14,FCGR1A,FCGR3A',
           # NK = 'NCAM1',
           GRANULOCYTE = 'CEACAM8',
           # 
           # #DC
           DCC = 'CD1C,CD83,CD209,THBD',
           DCP ='CLEC4C,NRP1,IL3RA',
           # 
           # #B CELL
           B = 'CD19',
           # B1 = 'MS4A1,CD27,SPN',
           # BNAIVE = 'CD22',
           # BMZ = '',
           BMATURE = 'IGHM',
           BPLASMA = 'IGHG1,IGHG2,IGHG3,SDC1',
           # BMEM = 'CD27,BCL2,NT5E',
           BREG = 'IL10')

for(i in h2){
  stringr::str_split(i,',')
}
h2 <- purrr::map(h2,function(i){
  strsplit(i,',')
})
saveRDS(h2,'searchF3.rds')


# cluster ------------------------------------------------------------------

for(i in names(searchF2)){
  print(i)
  print(searchF2[[i]][[1]][!(searchF2[[i]][[1]] %in% feats)])
}
