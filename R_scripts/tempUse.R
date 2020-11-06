#! /home/skelviper/anaconda3/envs/R/bin/R

# scripts for generate differnet types of pairs file : e.g.  interChromosome intraChromosome longRange interactions...

# change above to your own R path
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0 | args[1]=="-h") {
  stop("usage: Rscript processPairs.R input.pairs output.pairs", call.=FALSE)
}

library(tidyverse)

write_tsv(x=read_table2(args[1])%>%filter(X1==X4&abs(pos1-pos2)>20000),path=args[2])