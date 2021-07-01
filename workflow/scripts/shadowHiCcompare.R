#!/usr/bin/env Rscript

# Set .libPaths to last path (the Conda environemnt) to avoid conflicts
paths = .libPaths()
lastPath = paths[length(paths)]
.libPaths(lastPath)

library(HiCcompare)

args = commandArgs(trailingOnly=TRUE)
matrix1 = args[1]
matrix2 = args[2]
out = args[3]
chr = args[4]

data.table <- create.hic.table(read.table(matrix1), read.table(matrix2), chr=chr)
hic.table <- hic_loess(data.table, Plot=FALSE, Plot.smooth=FALSE)
hic.table = as.data.frame(hic.table)[, c('start1', 'start2', 'adj.M')]
write.table(hic.table, out, sep=' ', quote=FALSE, row.names=FALSE, col.names=FALSE)