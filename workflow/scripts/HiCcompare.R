#!/usr/bin/env Rscript

library(reshape2)
library(HiCcompare)
#pdf(NULL)

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


missingBins <- function(hic.table) {
  # Some intervals may be completely missing fomr bin1 or bin2 due to missing data.
  # Dummy values must be added here to ensure a square matrix is produced.
  missingBins = unique(c(
    setdiff(unique(hic.table$bin1), unique(hic.table$bin2)), 
    setdiff(unique(hic.table$bin1), unique(hic.table$bin2))))
  for (bin in missingBins) {
    hic.table[nrow(hic.table)+1, c("bin1","bin2")]  = list(bin, bin)
  }
  return(hic.table)
}


writeMatrix <- function(hic.table, out) {
  hic.table$bin1 = paste(hic.table$chr1, hic.table$start1, sep='-')
  hic.table$bin2 = paste(hic.table$chr2, hic.table$start2, sep='-')
  hic.table = missingBins(hic.table)
  homer <- dcast(hic.table, bin1 ~ bin2, value.var = "Z", fill = 0)
  rows <- homer[,1]
  homer[,1] <- NULL
  homer[lower.tri(homer)] <- t(homer)[lower.tri(homer)]
  homer = cbind.data.frame('HiCMatrix (directory=.)' = rows, 'Regions' = rows, homer)
  write.table(x = homer, file = out, quote = FALSE, sep = "\t", row.names = FALSE)
}


args = commandArgs(trailingOnly=TRUE)

outdir = args[1]
chr = args[2]
binsize = as.integer(args[3])
matrices = tail(args, -3)

numChanges = as.integer(2.539842873*10^-8 * binsize^2 - 2.871604938*10^-2 * binsize + 3617.620651)

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

    filter_params(hic.table, numChanges = numChanges, Plot = TRUE)
    png(filter_plot)

    hic.table <- hic_compare(hic.table, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE)
    png(compare_plot)

    hic.table$abs.adj.M = abs(hic.table$adj.M)
    hic.table$score = (hic.table$abs.adj.M / max(abs(hic.table$adj.M))) * 1000

    # Use absolute M for pyGenomeTracks, add normal M to last column to be read by links2interval script
    write.table(
      hic.table[,c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'abs.adj.M', 'score', 'adj.M', 'p.adj')],
      out_links, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
    
    
    writeMatrix(as.data.frame(hic.table), out_matrix)
  }
}





