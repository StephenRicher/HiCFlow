#!/usr/bin/env Rscript

# Set .libPaths to last path (the Conda environemnt) to avoid conflicts
paths = .libPaths()
lastPath = paths[length(paths)]
.libPaths(lastPath)

library(reshape2)
library(HiCcompare)

# Prevent scientific notation
options(scipen=999)

get_group <- function(path) {
  splitpath = gsub('\\..*$', '', basename(path))
  split = strsplit(splitpath, '-')[[1]]
  sample = split[1]
  return(sample)
}

args = commandArgs(trailingOnly=TRUE)

outdir = args[1]
qcdir = args[2]
chr = args[3]
start = as.integer(args[4])
end = as.integer(args[5])
binsize = as.integer(args[6])
suffix = args[7]
matrix1 = args[8]
matrix2 = args[9]

dir.create(qcdir, recursive = TRUE)

group1 = get_group(matrix1)
group2 = get_group(matrix2)

loess_plot = paste(qcdir, '/', group1, '-vs-', group2, '-loess-', suffix, '.png', sep = '')
outIF1 = paste(outdir, '/', group1, '-vs-', group2, '-adjIF1-', suffix, '.2d.txt', sep = '')
outIF2 = paste(outdir, '/', group1, '-vs-', group2, '-adjIF2-', suffix, '.2d.txt', sep = '')
outIF1SUTM = paste(outdir, '/', group1, '-vs-', group2, '-adjIF1-', suffix, '.sutm', sep = '')
outIF2SUTM = paste(outdir, '/', group1, '-vs-', group2, '-adjIF2-', suffix, '.sutm', sep = '')

### LOESS  ###
data.table <- create.hic.table(
  read.table(matrix1), read.table(matrix2), chr = chr, include.zeros = TRUE)

png(loess_plot)
hic.table <- hic_loess(data.table, Plot=TRUE, Plot.smooth=FALSE)
dev.off()

hic.table = as.data.frame(hic.table)

# Write adjusted IF1 to SUTM
write.table(
  hic.table[, c('start1', 'start2', 'adj.IF1')],
  outIF1SUTM, sep=' ', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Write adjusted IF2 to SUTM
write.table(
  hic.table[, c('start1', 'start2', 'adj.IF2')],
  outIF2SUTM, sep=' ', quote=FALSE, row.names=FALSE, col.names=FALSE)


# Create copy of data frame for lower triangle
hic.table2 = data.frame(hic.table)
names(hic.table2) = c("chr2", "start2",  "end2", "chr1", "start1", "end1", "IF1", "IF2", "D", "M", "adj.IF1", "adj.IF2", "adj.M",  "mc", "A")
# Remove the centre diagonal to avoid two copies
hic.table2 = hic.table2[(hic.table2$start2 != hic.table2$start1),]
# Combine upper and lower values
hic.table = rbind(hic.table, hic.table2)

hic.table['start1'] = hic.table['start1'] * binsize
hic.table['end1'] = hic.table['end1'] * binsize
hic.table['start2'] = hic.table['start2'] * binsize
hic.table['end2'] = hic.table['end2'] * binsize

# Write adjusted IF1 values
write.table(
  hic.table[,c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'adj.IF1')],
  outIF1, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Write adjusted IF2 values
write.table(
  hic.table[,c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'adj.IF2')],
  outIF2, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
