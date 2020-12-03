#!/usr/bin/env Rscript

library(hicrep)

get_sample <- function(path) {
  splitpath = gsub('\\..*$', '', basename(path))
  split = strsplit(splitpath, '-')[[1]]
  sample = paste(split[1], split[2], sep = '-')
  return(sample)
}

args = commandArgs(trailingOnly=TRUE)

outfile = args[1]
bin = as.integer(args[2])
start = as.integer(args[3])
end = as.integer(args[4])
matrix1 = args[5]
matrix2 = args[6]

h_hat = NA
scc = NA

if (file.info(matrix1)$size == 0) {
  print(paste(matrix1 ,"is empty - SCC not computed."))
} else if (file.info(matrix2)$size == 0) {
  print(paste(matrix2 ,"is empty - SCC not computed."))
} else {
  # Set max interaction as half capture region size, or 1,000,000bp
  max_interaction = min(1000000, as.integer((end - start)/2))
  
  m1 = read.table(matrix1)
  m2 = read.table(matrix2)
  
  # Calculate optimal smoothing parameter for a given region and bin
  h_hat = tryCatch({
    htrain(m1, m2, resol = bin, max = max_interaction, range = 0:20)
  }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  
  if (is.null(h_hat)) {
    print("Unable to calculate smoothing parameter - SCC note computed")
    h_hat = NA
    scc = NA
  } else {
    # Downsample matrices to the same sequencing depth
    min_sample = as.integer(min(sum(m1[,-c(1:3)]), sum(m2[,-c(1:3)])))
    
    m1_ds <- depth.adj(m1, min_sample, bin, out = 0)
    m2_ds <- depth.adj(m2, min_sample, bin, out = 0)
    
    # Format HiC matrix pairs and smooth
    pre_hic <- prep(m1_ds, m2_ds, resol = bin, h = h_hat, max = max_interaction)
    
    # Calculate SCC
    scc = get.scc(pre_hic, resol = bin, max = max_interaction)$scc
  }
}

hicrepResult <- data.frame(
  "sample1" = get_sample(matrix1), "sample2" = get_sample(matrix2),
  "h" = h_hat, "scc" = scc, stringsAsFactors=FALSE)
write.table(hicrepResult, outfile, sep=',', row.names=FALSE)