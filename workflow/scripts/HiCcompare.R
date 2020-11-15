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
qcdir = args[2]
chr = args[3]
start = as.integer(args[4])
end = as.integer(args[5])
binsize = as.integer(args[6])
fdr = as.double(args[7])
matrix1 = args[8]
matrix2 = args[9]

dir.create(qcdir)

group1 = get_group(matrix1)
group2 = get_group(matrix2)

loess_plot = paste(qcdir, '/', group1, '-vs-', group2, '-loess.png', sep = '')
filter_plot = paste(qcdir, '/', group1, '-vs-', group2, '-filter_params.png', sep = '')
compare_plot = paste(qcdir, '/', group1, '-vs-', group2, '-hicCompare.png', sep = '')
out_links = paste(outdir, '/', group1, '-vs-', group2, '.links', sep = '')
out_matrix = paste(outdir, '/', group1, '-vs-', group2, '.homer', sep = '')
out_sig = paste(outdir, '/', group1, '-vs-', group2, '-sig.homer', sep = '')
out_fdr = paste(outdir, '/', group1, '-vs-', group2, '-fdr.homer', sep = '')
out_medianAdjM = paste(outdir, '/', group1, '-vs-', group2, '-absZ.bedgraph', sep = '')

data.table <- create.hic.table(read.table(matrix1), read.table(matrix2), chr = chr)

png(loess_plot)
hic.table <- hic_loess(data.table, Plot = TRUE, Plot.smooth = FALSE)
dev.off()

# Number of changes is 1% of unique bins or 300, whichever higher
changes = as.integer(max(300, numBins(start, end, binsize) * 0.01))

# Very sparse matrices can trigger exceptions in filter params
png(filter_plot)
try(filter_params(hic.table, numChanges = changes, Plot = TRUE))
dev.off()

png(compare_plot)
hic.table <- hic_compare(hic.table, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE)
dev.off()

hic.table = as.data.frame(hic.table)

# Write mean adjM per interval as bedgraph
require(data.table)
medianAdjM = data.table(hic.table)
medianAdjM = medianAdjM[order(medianAdjM$start1), list(abs.Z = mean(abs(Z), na.rm = TRUE), adjM.count=length(adj.M)),
                        by=c("chr1","start1", "end1")]

# Ensure start and end bedgraph intervals extend to to full interval
if (min(medianAdjM$start1) > start) {
  medianAdjM = rbindlist(list(
    list(chr, start, min(medianAdjM$start1), 0, 0), medianAdjM))
}
if (max(medianAdjM$end1) < end) {
  medianAdjM = rbindlist(list(
    medianAdjM, list(chr, max(medianAdjM$end1), end, 0, 0)))
}
# Set all median values computed by less than 'n' values to 0
medianAdjM[medianAdjM$adjM.count < 10, 'abs.Z'] = 0
write.table(medianAdjM, out_medianAdjM, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

# Write matrix of logFC values
writeMatrix(hic.table, out_matrix, chr, start, end, binsize, 'adj.M')

# Write matrix of of significant FDR values
writeMatrix(hic.table[hic.table$p.adj <= fdr,], out_sig, chr, start, end, binsize, 'adj.M')

# Write matrix of signed inverted adjusted P values
hic.table$scale.p.adj = 1 - as.double(hic.table$p.adj)
hic.table$sign.p.adj = ifelse(hic.table$adj.M > 0, hic.table$scale.p.adj, -hic.table$scale.p.adj)

writeMatrix(hic.table, out_fdr, chr, start, end, binsize, "sign.p.adj")


hic.table$abs.adj.M = abs(hic.table$adj.M)
hic.table$score = (hic.table$abs.adj.M / max(abs(hic.table$adj.M))) * 1000

# Use absolute M for pyGenomeTracks, add normal M to last column to be read by links2interval script
write.table(
  hic.table[,c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'abs.adj.M', 'score', 'adj.M', 'p.adj')],
  out_links, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')


