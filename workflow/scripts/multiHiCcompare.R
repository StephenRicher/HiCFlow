#!/usr/bin/env Rscript

library(reshape2)
library(multiHiCcompare)
library(GenomicRanges)
# Prevent scientific notation
options(scipen=999)

split_name <- function(path) {
  splitpath = gsub('\\..*$', '', basename(path))
  split = strsplit(splitpath, '-')[[1]]
  return(list('group' = split[1], 'rep' = split[2]))
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

args = commandArgs(trailingOnly=TRUE)

outdir = args[1]
qcdir = args[2]
chr = args[3]
start = as.integer(args[4])
end = as.integer(args[5])
binsize = as.integer(args[6])
fdr = as.double(args[7])
matrices = args[-c(1:7)]

dir.create(qcdir)

args = vector("list", length(matrices))
groups = c()
n = 1
for (matrix in matrices) {
  name = split_name(matrix)
  groups = c(groups, name$group)
  # Substitute chr name with 99 as multiHiCcompare requires numeric.
  data = cbind(99, read.table(matrix))
  args[[n]] = data
  n = n + 1
}

hicexp = make_hicexp(data_list = args, groups = groups, zero.p = 0.8, 
                     A.min = 5, filter = TRUE, remove.regions = NULL)
hicexp <- cyclic_loess(hicexp, verbose = FALSE, parallel = FALSE, span = 0.2)

md_plot = paste(qcdir, '/', unique(groups)[1], '-vs-', unique(groups)[2], '-loess.png', sep = '')
png(md_plot)
MD_hicexp(hicexp)
dev.off()

hicexp <- hic_exactTest(hicexp, p.method = 'fdr', parallel = FALSE)

# Write matrix of logFC values
hic.table = as.data.frame(
  topDirs(hicexp, logfc_cutoff = 0, logcpm_cutoff = 0, p.adj_cutoff = 1, 
          D_cutoff = 0, alpha = 1, return_df = "pairedbed"))
out_matrix = paste(outdir, '/', unique(groups)[1], '-vs-', unique(groups)[2], '.homer', sep = '')
writeMatrix(hic.table, out_matrix, chr, start, end, binsize, "logFC")

# Write matrix of of significant FDR values
hic.table.sig = as.data.frame(
  topDirs(hicexp, p.adj_cutoff = fdr, return_df = 'pairedbed'))
out_sig = paste(outdir, '/', unique(groups)[1], '-vs-', unique(groups)[2], '-sig.homer', sep = '')
writeMatrix(hic.table.sig, out_sig, chr, start, end, binsize, "logFC")

# Write matrix of signed inverted adjusted P values
hic.table$scale.p.adj = 1 - as.double(hic.table$p.adj)
hic.table$sign.p.adj = ifelse(hic.table$logFC > 0, hic.table$scale.p.adj, -hic.table$scale.p.adj)
out_fdr = paste(outdir, '/', unique(groups)[1], '-vs-', unique(groups)[2], '-fdr.homer', sep = '')
writeMatrix(hic.table, out_fdr, chr, start, end, binsize, "sign.p.adj")


td <- topDirs(hicexp, logfc_cutoff = 1, logcpm_cutoff = 0.5, 
              p.adj_cutoff = fdr, return_df = 'pairedbed')
td$abs.logFC = abs(td$logFC)
td$score = (td$abs.logFC / max(td$abs.logFC)) * 1000
# Use absolute M for pyGenomeTracks, add normal M to last column to be read by links2interval script
out_links = paste(outdir, '/', unique(groups)[1], '-vs-', unique(groups)[2], '.links', sep = '')
write.table(
  td[,c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'abs.logFC', 'score', 'logFC', 'p.adj')],
  out_links, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')