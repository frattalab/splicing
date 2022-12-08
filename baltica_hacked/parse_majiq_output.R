#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(readr)
  library(dplyr)
  library(optparse)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "majiq/voila/*_voila.tsv",
    help = "Path with glob character to Majiq result files. [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "majiq/majiq_junctions.csv",
    help = "Path to output file [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cutoff"),
    default = 1.1,
    type = "double",
    help = "Discard junction with probability_non_changing < than
    --cutoff [default %default, i.e. no filter]"
  )
)
# enabling both baltica or snakemake input
if (exists("snakemake")) {
  opt <- list(
    input = snakemake@input,
    output = snakemake@output[[1]],
    cutoff = snakemake@params[["cutoff"]]
  )
  files <- opt$input
  file_names <- gsub(
    x = opt$input,
    pattern = "majiq/voila/(.+)_voila.tsv",
    replacement = "\\1"
  )
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
  files <- Sys.glob(opt$input)
  file_names <- gsub(
    x = files,
    replacement = "\\1",
    pattern = sub(x = opt$input, pattern = "\\*", replacement = "(.*)")
  )
}

# rename column names from majiq result due to the presence of spaces
read_majiq_out <- function(x) {
  read.table(
    x,
    sep = "\t",
    header = TRUE,
    col.names = c(
      "gene_name",
      "gene_id",
      "lsv_id",
      "mean_dpsi_per_lsv_junction",
      "probability_changing",
      "probability_non_changing",
      "ref_mean_psi",
      "alt_mean_psi",
      "lsv_type",
      "num_junctions",
      "num_exons",
      "de_novo_junctions",
      "chr",
      "strand",
      "junctions_coords",
      "exons_coords",
      "ir_coords",
      "ucsc_lsv_link"
    ),
    colClasses = c(
      "character",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character",
      "integer",
      "integer",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character"
    )
  )
}

res <- lapply(files, read_majiq_out)

names(res) <- file_names

res <- bind_rows(res, .id = "comparison") %>%
  separate_rows(
    "junctions_coords",
    "mean_dpsi_per_lsv_junction",
    "probability_changing",
    "probability_non_changing",
    "ref_mean_psi",
    "alt_mean_psi",
    sep = ";",
    convert = T
  )


junction_pattern <- "(\\d+)-(\\d+)"
junctions_coords <- str_match(
  res$junctions_coords, junction_pattern
)[, c(2, 3)]

res["start"] <- junctions_coords[, 1]
res["end"] <- junctions_coords[, 2]

message("Number of junctions output by Majiq ", nrow(res))
res <- dplyr::select(res, c(
  chr,
  start,
  end,
  strand,
  comparison,
  probability_changing,
  probability_non_changing,
  lsv_id,
  ref_mean_psi,
  alt_mean_psi
)) %>%
  filter(probability_non_changing < opt$cutoff) %>%
  mutate(method = "majiq") %>%
  arrange(probability_non_changing) %>%
  distinct(comparison, chr, start, end, strand, .keep_all = TRUE)


message("Number of junctions after filtering ", nrow(res))

write_csv(res, opt$output)