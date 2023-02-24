#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(readr)
  library(dplyr)
  library(optparse)
  library(GenomicRanges)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "leafcutter/*/*_cluster_significance.txt",
    help = "Path with glob character to Leafcutter result
    files. [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "leafcutter/leafcutter_junctions.csv",
    help = "Path to output file [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cutoff"),
    type = "double",
    default = 1.1,
    help = "Discard junction with p.adjust < than --cutoff
    [default %default, i.e. no filter]"
  )
)
# enabling both baltica or snakemake input
if (exists("snakemake")) {
  opt <- list(
    input = snakemake@input,
    output = snakemake@output[[1]],
    cutoff = snakemake@params[["cutoff"]]
  )
  files <- as.character(opt$input)
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
  files <- Sys.glob(opt$input)
}

file_names <- str_split(files, "/", simplify = T)
file_names <- file_names[, ncol(file_names) - 1]
message("Comparison names: ", paste0(file_names, collapse = "\t"))

message("Loading the cluster significance files")
cluster_sig_file <- files

cluster_sig <- lapply(
  cluster_sig_file,
  read.table,
  sep = "\t",
  header = TRUE
)

names(cluster_sig) <- file_names
cluster_sig <- bind_rows(cluster_sig, .id = "comparison")

message("Loading the effect sizes files")

effec_size_files <- gsub(
  x = files,
  pattern = "cluster_significance",
  replacement = "effect_sizes"
)

es <- lapply(
  effec_size_files,
  read.table,
  sep = "\t",
  header = T,
  colClasses = c("character", "double", "double", "double", "double")
)

names(es) <- file_names
es <- bind_rows(es, .id = "comparison")

# parse the intron column for merging
es$intron <- str_replace_all(es$intron, "_", ":")

intron <- str_split(es$intron, pattern = ":", simplify = T)
colnames(intron) <- c(
  "chr", "start", "end", "clu", "clu_number", "strand"
)
intron <- as_tibble(intron)

es <- bind_cols(es, intron)
# cluster will be the pivot for merging
es$cluster <- as.character(
  str_glue_data(es, "{chr}:{clu}_{clu_number}_{strand}")
)
message("Merging tables")

res <- inner_join(es, cluster_sig, by = c("comparison", "cluster"))
res$chr <- gsub("chr", "", res$chr)

res <- select(res, -c("clu", "clu_number"))
# create a unique junction column for each row
message("Number of junctions output by Leafcutter ", nrow(res))
res <- res %>%
  filter(p.adjust < opt$cutoff) %>%
  mutate(method = "Leafcutter") %>%
  arrange(p.adjust) %>%
  distinct(comparison, chr, start, end, strand, .keep_all = TRUE)


message("Number of junctions after filtering ", nrow(res))
write_csv(res, opt$output)