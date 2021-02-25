# libraries ---------------------------------------------------------------

library(rvest)
library(tidyverse)

# parse17-20 -----------------------------------------------------------------

url1 <- "http://fcb.ycga.yale.edu:3010/cTNWq4WM7RGfToedjNA05zKEwQmIb/031620/"

links <- read_html(url1) %>% 
  html_nodes('a') %>% 
  html_attr('href')
links <- links[str_detect(links,'.*/$')][-1]

filetype <- 'filtered_contig_annotations.csv'
# url2 <- str_c(url1,links,filetype)

for (link in links) {
  data <- data.table::fread(paste0(url1,link,filetype))
  write.csv(data, file=paste0('data_tcr/',
    str_extract(link,'.*(?=-6_VHT_cellranger)'), '.csv'))
}

# parse21-24 ---------------------------------------------------------------

url1 <- "http://fcb.ycga.yale.edu:3010/GDQlMrcQ2UPFclPUH5Z9Pqp9_Nmmi/111620/"

links <- read_html(url1) %>% 
  html_nodes('a') %>% 
  html_attr('href')
links <- links[str_detect(links,'.*ranger/$')]

for (link in links) {
  data <- data.table::fread(paste0(url1,link,filetype))
  write.csv(data, file=paste0('data_tcr/',
                              str_extract(link,'.*(?=-6VDJ)'), '.csv'))
}

