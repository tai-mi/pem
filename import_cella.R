# library -----------------------------------------------------------------

library(tidyverse)
library(CellaRepertorium)


# Load ---------------------------------------------------------------

fList <- list.files('tcr_data')
tempData <- data.frame()

for (file in fList) {
  tempData <- rbind(tempData,
                    data.table::fread(paste0('tcr_data/',file))[,-1])
}


# QC ----------------------------------------------------------------------

print(paste(
  ifelse(all(h1$contig_tbl$high_confidence),
         'All high confidence,',
         paste0(100*sum(h1$contig_tbl$high_confidence)/length(h1$contig_tbl$high_confidence),'% of cells are high confidence,')),
  ifelse(all(h1$contig_tbl$is_cell),
         'all from a cell,',paste0(100*sum(h1$contig_tbl$is_cell)/length(h1$contig_tbl$is_cell),'% of contigs are from cells,')),
  ifelse(all(h1$contig_tbl$full_length),
         'all full length,',
         paste0(100*sum(h1$contig_tbl$full_length)/length(h1$contig_tbl$full_length),
                '% of contigs are full length,')),'and',
  ifelse(all(h1$contig_tbl$productive),'all productive',
         paste0(100*sum(h1$contig_tbl$productive)/length(h1$contig_tbl$productive),
                '% are productive'))
))
