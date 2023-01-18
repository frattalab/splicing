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
    default = 1.1,
    help = "Discard junction with FDR < than --cutoff [default %default]"
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

#' Process for alternative splice site rMATs output files
#' @param df dataframe from rMATs
#' @param start name of the start column
#' @param end name of the end column
#' @param type flag for splice junction type
#' @param FDR is the FDR cutoff
#' @return GenomicRange of the selected SJ
#' @export
#'
process_RMATS_ass <- function(df, start = "flankingES", end = "shortES", type, FDR = 0.05) {
  df <- df %>%
    dplyr::filter(FDR < !!FDR) %>%
    dplyr::select(
      chr, !!start, !!end, strand, comparison, FDR, IncLevelDifference
    ) %>%
    dplyr::rename(c(start = !!start, end = !!end)) %>%
    mutate(start = pmin(start, end), end = pmax(start, end))

  df$"type" <- type

  df
}


#' Process for exon skipping and intron retention r
#' MATs output files
#' @param df dataframe from rMATs
#' @param start name of the start column
#' @param end name of the end column
#' @param type flag for splice junction type
#' @param FDR is the FDR cutoff
#' @return GenomicRange of the selected SJ
#' @export
#'
process_RMATS <- function(df, start, end, type, FDR = 0.05) {
  df <- df %>%
    dplyr::filter(FDR < !!FDR) %>%
    dplyr::select(
      chr, !!start, !!end, strand, comparison,
      FDR, IncLevelDifference
    ) %>%
    dplyr::rename(c(start = !!start, end = !!end))

  df$type <- type

  df
}

split_path <- strsplit(files, "/", fixed = T)
comparison <- sapply(split_path, `[[`, length(split_path[[1]]) - 1)

res <- split(files, comparison)

get_rmats_coord <- function(.files, .comparison) {
  message("Processing files for ", .comparison)
  x <- lapply(.files, read.delim, "\t", header = TRUE)

  lapply(seq_along(.files), function(i) {
    message(.files[[i]], ": has ", nrow(x[[i]]), " entries")
  })
  .split_files <- strsplit(.files, "/", fixed = T)
  as_type <- gsub(
    x = sapply(.split_files, `[[`, length(.split_files[[1]])),
    pattern = ".MATS.JC.txt", replacement = ""
  )

  names(x) <- as_type
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

  res <- bind_rows(
    lst(es_ssj, es_isj, ir, a5ss_SSJ, a5ss_ISJ, a3ss_SSJ, a3ss_ISJ)
  )
  res$method <- "rmats"
  res
}

message("Loading and procesing rMATs files")
res <- lapply(setNames(names(res), names(res)), function(x) {
  get_rmats_coord(res[[x]], x)
})
res <- bind_rows(res)

# Force remove chr from seqnames
res <- res %>%
  mutate(
    chr = gsub(x = res$chr, replacement = "", pattern = "chr")
  ) %>%
  arrange(FDR) %>%
  distinct(comparison, chr, start, end, strand, .keep_all = TRUE)

message("Number of junctions after filtering ", nrow(res))
write_csv(res, opt$output)