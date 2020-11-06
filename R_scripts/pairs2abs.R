#! /home/skelviper/anaconda3/envs/R/bin/R

# change above to your own R path
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0 | args[1]=="-h") {
  stop("usage: Rscript processPairs.R chrLength input.pairs output.abs.pairs", call.=FALSE)
}

library(tidyverse)

chrLength <- read_table2(args[1],col_names = FALSE)
names(chrLength) <- c("chr","length","absolutePos")

#别数了，就是23条。
chrlist = c('chr14','chr5','chr15','chr10','chr18','chr7','chr6','chr9','chr16','chr19','chrX','chr8','chr4','chr12','chr17','chr2','chr11','chr20','chr3','chr1','chr21','chr13','chr22')

pairsfile <- read_table2(args[2]) %>% filter(X4 %in% chrlist)
names(pairsfile) <- c("readID","chr1","pos1","chr2","pos2")
pairsfile <- pairsfile  %>% group_by(chr1,chr2,pos1,pos2,readID)

options(warn=-1)
absPos <- pairsfile %>% mutate(pos1 = round(pos1 + filter(chrLength,chr==chr1)[[1,3]]),pos2 = round(pos2 + filter(chrLength,chr==chr2)[[1,3]]))

write_tsv(x=absPos,path=args[3])


