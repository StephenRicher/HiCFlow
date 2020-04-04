#!/usr/bin/env Rscript

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

args = commandArgs(trailingOnly=TRUE)

outdir = args[1]
chr = args[2]
binsize = as.integer(args[3])
matrices = tail(args, -3)

# Remove any files that may be empty
for (matrix in matrices) {
  if (empty_file(matrix)) {
    matrices = matrices[matrices != matrix]
  }
}

numChanges = as.integer(2.539842873*10^-8 * binsize^2 - 2.871604938*10^-2 * binsize + 3617.620651)

# Note this script reruns pairwise groups to more easily define snakemake output

for (matrix1 in matrices) {
  for (matrix2 in matrices) {
    
    group1 = get_group(matrix1)
    group2 = get_group(matrix2)
    
    loess_plot = paste(outdir, '/', group1, '-vs-', group2, '-loess.png', sep = '')
    filter_plot = paste(outdir, '/', group1, '-vs-', group2, '-filter_params.png', sep = '')
    compare_plot = paste(outdir, '/', group1, '-vs-', group2, '-hicCompare.png', sep = '')
    out_links = paste(outdir, '/', group1, '-vs-', group2, '.links', sep = '')
    
    system(paste('touch', loess_plot))
    system(paste('touch', compare_plot))
    system(paste('touch', out_links))

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
      hic.table[,c('chr1', 'start1', 'start2', 'chr2', 'start2', 'end2', 'abs.adj.M', 'score', 'adj.M', 'p.adj')],
      out_links, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
  }
}




