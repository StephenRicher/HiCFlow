#!/usr/bin/env Rscript

library(reshape2)
library(HiCcompare)


addMissingIntervals <- function(hic.table, start, end, binsize) {
  if (nrow(hic.table) == 0) {
    intervals = seq(start, end, binsize)
  } else {
    intervals = seq(min(hic.table$start1), max(hic.table$start1), binsize)
  }
  # Add a 2 bin buffer to each side
  intervals = sort(unique(c(seq(min(intervals), start - (2*binsize), -binsize),
                            intervals, seq(max(intervals), end + (2*binsize), binsize))))
  intervals = intervals[!intervals < 0]
  for (place in c('start1', 'start2')) {
    for (interval in setdiff(intervals, unique(hic.table[,place]))) {
      hic.table[nrow(hic.table) + 1, c("start1","start2")]  = list(interval, interval)
    }
  }
  return(hic.table)
}


writeMatrix <- function(hic.table, out, chr, start, end, binsize, var) {
  hic.table = addMissingIntervals(hic.table, start, end, binsize)
  homer <- dcast(hic.table, start1 ~ start2, value.var = var, fill = 0)
  homer[,1] = paste(chr, homer[,1], sep='-')
  rows = homer[,1]
  homer[,1] = NULL
  colnames(homer) = rows
  homer[lower.tri(homer)] <- t(homer)[lower.tri(homer)]
  homer = cbind.data.frame('HiCMatrix (directory=.)' = rows, 'Regions' = rows, homer)
  write.table(x = homer, file = out, quote = FALSE, sep = "\t", row.names = FALSE)
}


numBins <- function(start, end, binsize) {
  intervals = length(seq(start, end, binsize))
  total_bins = intervals^2
  lower_half = (total_bins - intervals)/2
  unique_bins =  total_bins - lower_half
  return(unique_bins)
}

args = commandArgs(trailingOnly=TRUE)

matrix1 = args[1]
matrix2 = args[2]
out = args[3]
chr = args[4]
start = as.integer(args[5])
end = as.integer(args[6])
binsize = as.integer(args[7])

data.table <- create.hic.table(read.table(matrix1), read.table(matrix2), chr=chr)

hic.table <- hic_loess(data.table, Plot=FALSE, Plot.smooth=FALSE)

hic.table = as.data.frame(hic.table)

# WRITE TO SUTM INSTEAD #

# Write matrix of logFC values
writeMatrix(hic.table, out, chr, start, end, binsize, 'adj.M')
