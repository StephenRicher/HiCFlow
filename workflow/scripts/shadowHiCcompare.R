#!/usr/bin/env Rscript

library(reshape2)
library(HiCcompare)
#pdf(NULL)

# Prevent scientific notation
options(scipen=999)

get_group <- function(path) {
  splitpath = gsub('\\..*$', '', basename(path))
  split = strsplit(splitpath, '-')[[1]]
  sample = split[1]
  return(sample)
}


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

outdir = args[1]
chr = args[2]
start = as.integer(args[3])
end = as.integer(args[4])
binsize = as.integer(args[5])
suffix = args[6]
matrix1 = args[7]
matrix2 = args[8]

print(matrix1)
print(matrix2)
group1 = get_group(matrix1)
group2 = get_group(matrix2)

out_matrix = paste(outdir, '/', group1, '-vs-', group2, '-', suffix, '.homer', sep = '')

data.table <- create.hic.table(read.table(matrix1), read.table(matrix2), chr = chr)

hic.table <- hic_loess(data.table, Plot = FALSE, Plot.smooth = FALSE)

# Number of changes is 1% of unique bins or 300, whichever higher
changes = as.integer(max(300, numBins(start, end, binsize) * 0.01))

# Very sparse matrices can trigger exceptions in filter params
try(filter_params(hic.table, numChanges = changes, Plot = FALSE))

hic.table <- hic_compare(hic.table, adjust.dist = TRUE, p.method = 'fdr', Plot = FALSE)

hic.table = as.data.frame(hic.table)

# Write matrix of logFC values
writeMatrix(hic.table, out_matrix, chr, start, end, binsize, 'adj.M')
