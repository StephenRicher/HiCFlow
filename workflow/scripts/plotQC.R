#!/usr/bin/env Rscript

library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(readr)


args = commandArgs(trailingOnly=TRUE)
dir = args[1]
if (!endsWith(dir, '/')) {
  dir = paste(dir, '/',  sep = '')
}
width = 10
height = (9/16) * width

# Merge input files into single CSV
data <- args[-1] %>%
  lapply(read_tsv) %>%
  bind_rows

data$replicate = gsub("^.*-", "", data$sample)
data$sample_group = gsub("-.*$","",data$sample)

trans_stats = table(pull(data[data$interaction_type == "trans", "orientation"]),
                    pull(data[data$interaction_type == "trans", "sample"]))
trans_stats = sweep(trans_stats,2,colSums(trans_stats),`/`)
write.csv(trans_stats,
          file = paste(dir,"trans_stats.csv", sep = ""))


insert_size = ggplot(data[data$interaction_type == "cis", ],
                     aes(x = insert_size, colour = orientation)) +
  # Normalize each facet panel seperately by total count of that facet group
  geom_freqpoly(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), bins = 100) +
  geom_vline(xintercept = 1000, linetype = "dotted", color = "#1B9E77") +
  facet_grid(sample_group ~ replicate, scales = "free_y") +
  scale_x_continuous(trans = 'log10') +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  labs(title = NULL, subtitle = NULL, tag = NULL,
       x = "Insert size (bp)", y = "Frequency", caption = NULL) +
  theme_pubr(legend = "bottom") +
  guides(colour = guide_legend(nrow = 1))
ggsave(filename = paste(dir,"insert_size_frequency.png", sep = ""),
       plot = insert_size,
       dpi = 300, width = width, height = height)
rm(insert_size)

ditag_length = ggplot(data[, c("ditag_length", "sample")],
                      aes(x = ditag_length, colour = sample)) +
  geom_density() +
  geom_vline(xintercept = 1000, linetype = "dotted", color = "#1B9E77") +
  geom_vline(xintercept = 100, linetype = "dotted", color = "#1B9E77") +
  scale_x_continuous(trans = 'log10', limits = c(NA, 10000)) +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  labs(title = NULL, subtitle = NULL, tag = NULL,
       x = "Ditag length (bp)", y = "Frequency",
       caption = NULL) +
  theme_pubr(legend = "bottom") +
  guides(colour = guide_legend(nrow = 1))
ggsave(filename = paste(dir,"ditag_length.png", sep = ""),
       plot = ditag_length,
       dpi = 300, width = width, height = height)
