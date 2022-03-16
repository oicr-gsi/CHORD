#!/usr/bin/env Rscript

library(mutSigExtractor) #is this necessary?
library(CHORD)
#library(data.table)
#library(tidyverse)

library(BSgenome.Hsapiens.UCSC.hg38) #reference

args = commandArgs(trailingOnly=TRUE)

snvFile_loc <- args[1]
structuralFile_loc <- args[2]
  
structuraldf <- read.table(structuralFile_loc)[,c(4:5)]
names(structuraldf) <- c("sv_type", "sv_len") 
#package documentation says col. names aren't necessary, but had trouble running without _these_ names


contexts <- extractSigsChord(
  vcf.snv = snvFile_loc,
  df.sv = structuraldf,
  ref.genome=BSgenome.Hsapiens.UCSC.hg38,
  verbose=T
)

chord_output <- chordPredict(contexts, do.bootstrap=T, verbose=T)

write.table(
  chord_output,
  file = paste(snvFile_loc,".CHORD.hrd.txt")
)
