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
    help = "Path with glob character to Leafcutter result files. [default %default]",
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
    default = 0.05,
    help = "Discard junction with p.adjust < than --cutoff [default %default]"
  )
)
# enabling both baltica or snakemake input
if (exists('snakemake')) {
    opt <- list(
      input = snakemake@input,
      output = snakemake@output[[1]],
      cutoff = snakemake@params[['cutoff']]
      )
    files <- opt$input

  } else {

    opt <- parse_args(OptionParser(option_list = option_list))
    files <- Sys.glob(opt$input)
  }

file_names = str_split(files, '/', simplify = T)
file_names = file_names[, ncol(file_names) - 1  ] 
message('Comparison names: ', paste0(file_names, collapse='\t') )

message("Loading the cluster significance files")
cluster_sig_file <- files # Sys.glob(file.path(out.path, '/*/*_cluster_significance.txt'))

cluster_sig <- lapply(
  cluster_sig_file,
  read_tsv,
  col_types = c(
    cluster = col_character(),
    status = col_character(),
    loglr = col_double(),
    df = col_double(),
    p = col_double(),
    p.adjust = col_double(),
    genes = col_character()
  )
)

names(cluster_sig) <- file_names
cluster_sig <- bind_rows(cluster_sig, .id = 'comparison')

message("Loading the effect sizes files")
effec_size_files <- gsub(x = files, pattern = 'cluster_significance', replacement = 'effect_sizes')#  Sys.glob(file.path(out.path, '/*/*effect_sizes.txt'))
es <- lapply(
  effec_size_files,
  read_tsv,
  col_names = c('intron', 'logef', 'ref_psi', 'alt_psi', 'deltapsi'),
  skip = 1,
  col_types = c(
    .default = col_double(),
    intron = col_character(),
    logef = col_double(),
    deltapsi = col_double()
  )
)
names(es) <- file_names
es <- bind_rows(es, .id = 'comparison')

# parse the intron column for merging
es$intron <- str_replace_all(es$intron, '_', ':')

# intron <- read_delim(es$intron, delim = ':', col_names = c('chr', 'start', 'end', 'clu', 'clu_number', 'strand'))
# es <- bind_cols(es, intron)
es <- es |> separate(intron, into = c('chr', 'start', 'end', 'clu', 'clu_number', 'strand'),sep =  ':')
# cluster will be the pivot for merging
es$cluster <- as.character(str_glue_data(es, "{chr}:{clu}_{clu_number}_{strand}"))
message("Merging tables")

res <- inner_join(es, cluster_sig, by = c('comparison', 'cluster'))
res$chr <- gsub('chr', '', res$chr)
# add ranks for psi per cluster
# res <- res %>%
#   arrange(comparison, cluster, ref_psi) %>%
#   group_by(comparison, cluster) %>%
#   mutate(is_canonical = row_number() == 1) %>%
#   ungroup()

res <- select(res, -c('clu', 'clu_number'))
# create a unique junction column for each row
message('Number of junctions output by Leafcutter ', nrow(res))

res <- res %>%
  filter(p.adjust < opt$cutoff) %>%
  mutate(method = 'Leafcutter')

message('Number of junctions after filtering ', nrow(res))
write_csv(res, opt$output)
