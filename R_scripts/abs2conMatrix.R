#! /home/skelviper/anaconda3/envs/R/bin/R

# change above to your own R path
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0 | args[1]=="-h") {
  stop("usage: Rscript processPairs.R resolution input.pairs output.abs.pairs", call.=FALSE)
}

library(tidyverse)

resolution <- as.numeric(args[1])
pairsfile <- read_table2(args[2])%>% select(1:5)

contact_matrix <- pairsfile %>% mutate(pos1=round(pos1 / resolution),pos2=round(pos2/resolution)) %>% select(pos1,pos2) %>% group_by(pos1,pos2) %>% summarise(count = n())

write_tsv(x=contact_matrix,path=args[3],col_names=FALSE)

