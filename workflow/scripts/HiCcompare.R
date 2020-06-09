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

empty_file <- function(file) {
  if (!(file.exists(file)) || (file.info(file)$size == 0)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


addMissingIntervals <- function(hic.table, start, end, binsize) {
  intervals = seq(min(hic.table$start1), max(hic.table$start1), binsize)
  intervals = sort(unique(c(seq(min(intervals), start - binsize, -binsize), 
                            intervals, seq(max(intervals), end + binsize, binsize))))
  for (place in c('start1', 'start2')) {
    for (interval in setdiff(intervals, unique(hic.table[,place]))) {
      hic.table[nrow(hic.table) + 1, c("start1","start2")]  = list(interval, interval)
    }
  }
  return(hic.table)
}


writeMatrix <- function(hic.table, out, chr, start, end, binsize) {
  hic.table = addMissingIntervals(hic.table, start, end, binsize)
  homer <- dcast(hic.table, start1 ~ start2, value.var = "Z", fill = 0)
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
matrices = tail(args, -5)

#numChanges = as.integer(2.539842873*10^-8 * binsize^2 - 2.871604938*10^-2 * binsize + 3617.620651)

# Note this script reruns pairwise groups to more easily define snakemake output
for (matrix1 in matrices) {
  for (matrix2 in matrices) {

    group1 = get_group(matrix1)
    group2 = get_group(matrix2)

    loess_plot = paste(outdir, '/', group1, '-vs-', group2, '-loess.png', sep = '')
    filter_plot = paste(outdir, '/', group1, '-vs-', group2, '-filter_params.png', sep = '')
    compare_plot = paste(outdir, '/', group1, '-vs-', group2, '-hicCompare.png', sep = '')
    out_matrix = paste(outdir, '/', group1, '-vs-', group2, '.homer', sep = '')
    out_links = paste(outdir, '/', group1, '-vs-', group2, '.links', sep = '')

    for (file in c(loess_plot, filter_plot, compare_plot, out_links, out_matrix)) {
        system(paste('touch', file))
    }

    data.table <- create.hic.table(read.table(matrix1), read.table(matrix2), chr = chr)

    hic.table <- hic_loess(data.table, Plot = TRUE, Plot.smooth = FALSE)
    png(loess_plot)
    dev.off()
    
    # Number of changes is 1% of unique bins or 300, whichever higher
    changes = as.integer(max(300, numBins(start, end, binsize) * 0.01))
    
    # Very sparse matrices can trigger exceptions in filter params
    err = tryCatch(
      expr = {
        filter_params(hic.table, numChanges = changes, Plot = TRUE)
      }, 
      error = function(err) {
        return(NULL)
      }
    )
    if (is.null(err)) {
      print(paste('Skipping', matrix1, matrix2))
      next
    }
    
    png(filter_plot)
    dev.off()
    hic.table <- hic_compare(hic.table, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE)
    png(compare_plot)
    dev.off()
    hic.table$abs.adj.M = abs(hic.table$adj.M)
    hic.table$score = (hic.table$abs.adj.M / max(abs(hic.table$adj.M))) * 1000

    # Use absolute M for pyGenomeTracks, add normal M to last column to be read by links2interval script
    write.table(
      hic.table[,c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'abs.adj.M', 'score', 'adj.M', 'p.adj')],
      out_links, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
    
    hic.table = as.data.frame(hic.table)
    writeMatrix(hic.table, out_matrix, chr, start, end, binsize)
  }
}





