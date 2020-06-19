#!/usr/bin/env Rscript

library(hicrep)
library(pheatmap)
library(ggplot2)
library(viridisLite)

get_sample <- function(path) {
  splitpath = gsub('\\..*$', '', basename(path))
  split = strsplit(splitpath, '-')[[1]]
  sample = paste(split[1], split[2], sep = '-')
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

outfile = args[1]
bin = as.integer(args[2])
start = as.integer(args[3])
end = as.integer(args[4])
matrices = tail(args, -4)

# Remove any files that may be empty
for (matrix in matrices) {
  if (empty_file(matrix)) {
    matrices = matrices[matrices != matrix]
  }
}

# Initialise dataframe to store values
scc_values <- data.frame("sample1" = character(0), "sample2" = character(0),
                         "h" = numeric(0), "scc" = numeric(0),
                         stringsAsFactors=FALSE)

# Set max interaction as half capture region size, or 1,000,000bp
max_interaction = min(1000000, as.integer((end - start)/2))

h_hat = NULL
sample_combinations = c()
for (matrix1 in matrices) {
  m1 = read.table(matrix1)
  for (matrix2 in matrices) {
    if (matrix1 != matrix2) {
      # Don't rerun equivalent combinations
      combo = paste(sort(c(matrix1, matrix2)), collapse = ":")
      if (combo %in% sample_combinations) {
        next
      }
      sample_combinations = c(sample_combinations, combo)
      
      m2 = read.table(matrix2)
      # Calculate optimal smoothing parameter for a given region and bin
      h_hat = tryCatch({
        htrain(m1, m2, resol = bin, max = max_interaction, range = 0:20)
      }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
      
      if (is.null(h_hat)) {
        print(paste("Skipping", matrix1, matrix2))
        scc_values[nrow(scc_values) + 1,] = list(
          get_sample(matrix1), get_sample(matrix2), NA, NA)
        next
      }

      # Downsample matrices to the same sequencing depth
      min_sample = as.integer(min(sum(m1[,-c(1:3)]), sum(m2[,-c(1:3)])))
      
      m1_ds <- depth.adj(m1, min_sample, bin, out = 0)
      m2_ds <- depth.adj(m2, min_sample, bin, out = 0)

      # Format HiC matrix pairs and smooth
      pre_hic <- prep(m1_ds, m2_ds, resol = bin, h = h_hat, max = max_interaction)

      # Calculate SCC
      SCC.out = get.scc(pre_hic, resol = bin, max = max_interaction)

      scc_values[nrow(scc_values) + 1,] = list(
        get_sample(matrix1), get_sample(matrix2), h_hat, SCC.out$scc)

    } else {
      scc_values[nrow(scc_values) + 1,] = list(
        get_sample(matrix1), get_sample(matrix2), NA, NA)
    }

  }
}

if (nrow(scc_values) == 0) {
  quit()
}

vals = sort(unique(c(as.character(scc_values$sample1), as.character(scc_values$sample2))))
nm = matrix(NA, nrow=length(vals), ncol=length(vals), dimnames=list(vals, vals))
nm[as.matrix(scc_values[, c('sample1', 'sample2')])] <- scc_values[, 'scc']
nm[as.matrix(scc_values[, c('sample2', 'sample1')])] <- scc_values[, 'scc']

# Check if any non-NA values are less than 0
if (any(nm[!is.na(nm)] < 0)) {
  print('Negative correlations detected - this is usually due to very low sparsity and so will be excluded.')
  nm[nm < 0] = NA
}

# Remove any rows or columns where all are NA - pheatmap won't plot this
nm = nm[rowSums(is.na(nm)) != ncol(nm), colSums(is.na(nm)) != nrow(nm)]

heatmap = pheatmap(nm,
                   color = viridis(100), display_numbers = TRUE,
                   number_color = 'red',
                   fontsize = 8, fontsize_number = 10,
                   angle_col = 45, treeheight_col = 0)
ggsave(filename = outfile, plot = heatmap, width = 8, height = 5)

