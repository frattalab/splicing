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
    default = "rmats/*/*.MATS.JC.txt",
    help = "Path with glob character to Leafcutter result files.
    [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "rmats/rmats_junctions.csv",
    help = "Path to output file [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cutoff"),
    type = "double",
    default = 0.05,
    help = "Discard junction with FDR < than --cutoff [default %default]"
  )
)
# enabling both baltica or snakemake input
if (exists("snakemake")) {
  snakemake@source("utils.R")

  opt <- list(
    input = snakemake@input,
    output = snakemake@output[[1]],
    cutoff = snakemake@params[["cutoff"]]
  )
  files <- opt$input
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
  files <- Sys.glob(opt$input)
  source("~/Baltica/scripts/utils.R")
}

res <- split(files, str_split(files, "/", simplify = T)[, 2])

get_rmats_coord <- function(.files, .comparison) {
  message("Processing files for ", .comparison)
  x <- suppressWarnings(
    lapply(.files, readr::read_delim, "\t", col_types = c(
      .default = col_double(),
      GeneID = col_character(),
      geneSymbol = col_character(),
      chr = col_character(),
      strand = col_character(),
      IJC_SAMPLE_1 = col_character(),
      SJC_SAMPLE_1 = col_character(),
      IJC_SAMPLE_2 = col_character(),
      SJC_SAMPLE_2 = col_character(),
      IncLevel1 = col_character(),
      IncLevel2 = col_character()
    ))
  )
  names(x) <- str_split(.files, pattern = "[/.]", simplify = T)[, 3]
  for (i in names(x)) {
    x[[i]]$comparison <- .comparison
  }

  es_ssj <- process_RMATS(
    x$SE, "upstreamEE", "downstreamES", "ES_SJ",
    FDR = opt$cutoff
  )
  es_isj <- process_RMATS(
    x$SE, "upstreamEE", "exonStart_0base", "EI_SJ",
    FDR = opt$cutoff
  )
  ir <- process_RMATS(
    x$RI, "upstreamEE", "downstreamES", "IR",
    FDR = opt$cutoff
  )
  a5ss_SSJ <- process_RMATS_ass(
    x$A5SS, "longExonEnd", "flankingES", "A5SS_SSJ",
    FDR = opt$cutoff
  )
  a5ss_ISJ <- process_RMATS_ass(
    x$A5SS, "shortEE", "flankingES", "A5SS_ISJ",
    FDR = opt$cutoff
  )
  a3ss_SSJ <- process_RMATS_ass(
    x$A3SS, "flankingEE", "longExonStart_0base", "A5SS_SSJ",
    FDR = opt$cutoff
  )
  a3ss_ISJ <- process_RMATS_ass(
    x$A3SS, "flankingES", "shortES", "A5SS_ISJ",
    FDR = opt$cutoff
  )

  gr <- c(es_ssj, es_isj, ir, a5ss_SSJ, a5ss_ISJ, a3ss_SSJ, a3ss_ISJ)
  gr$method <- "rmats"

  gr
}

message("Loading and procesing rMATs files")
res <- lapply(setNames(names(res), names(res)), function(x) {
  get_rmats_coord(res[[x]], x)
})
res <- as(res, "GRangesList")
res <- unlist(res)

message("Number of junctions after filtering ", length(res))
write_csv(data.frame(res), opt$output)