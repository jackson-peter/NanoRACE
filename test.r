library(tidyverse)
library(data.table)

flist <- list.files(path="/home/jpeter/DATA/FLEPseq/JP_TestRACE/3_PolyA/", pattern = "*.read_info.result.merged.txt", full.names = T)
names(flist)=lapply(flist, basename)

tot <- rbindlist(lapply(flist, fread), idcol = "barcode")
tot <- tot%>% mutate(barcode=gsub(".nanopore.read_info.result.merged.txt", "", barcode))

tot1 <- tot%>% group_by(barcode, type) %>%
  summarise(count_type=n())

ggplot(tot1, aes(x=type, y=count(type))) + geom_col() + facet_wrap(~barcode)
