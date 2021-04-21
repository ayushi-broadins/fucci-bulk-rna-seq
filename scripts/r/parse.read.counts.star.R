#!/usr/bin/env Rscript

parse.star.counts <- function(fname, strand){
  if(strand == "both"){
    count.col <- "V2"
  } else if(strand == "forward"){
    count.col <- "V3"
  } else if(strand == "reverse"){
    count.col <- "V4"
  }
  
  x <- fread(fname)
  x <- x[,c("V1",count.col), with=F] # keep only counts with specified strandedness
  sample <- strsplit(fname, "/")[[1]]
  sample <- sample[length(sample)-1]
  setnames(x, c("V1", count.col), c("gene_id", sample))
  return(x)
}


library(data.table)

args <- commandArgs(trailingOnly=T)

infiles <- strsplit(args[1], ";")[[1]]
strand <- args[2]
outfile <- args[3]

counts.list <- lapply(infiles, parse.star.counts, strand = strand)
counts <- Reduce(function(x,y) merge(x,y,by="gene_id"), counts.list)

write.table(counts, outfile, sep="\t", col.names=T, row.names=F, quote=F)
