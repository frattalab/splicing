#!/usr/bin/env Rscript
# Title     : annotate_SJ from many DJU methods
# Objective : To annotate SJ called DS by many DJU methods with gene
# and transcript level information
# Created by: Thiago Britto-Borges (thiago.brittoborges@uni-heidelberg.de)
# Created on: 30.04.20
suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(rtracklayer)
  library(openxlsx)
  library(readr)
  library(dplyr)
  library(stringr)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Path to parsed DJU result", ,
    metavar = "character",
    default = "results/SJ_annotated_assigned.csv"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    help = "Path to the output file",
    metavar = "character",
    default = "results/SJ_annotated_assigned_simple.xlsx"
  )
)

if (exists("snakemake")) {
  opt <- list(
    input = snakemake@input[[1]],
    output = snakemake@output[[1]]
  )
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
}

df <- read.csv(opt$input)

simplify <- function(x, remove = c()) {
  x %>%
    str_split(";") %>%
    lapply(., function(x) {
      setdiff(x, remove)
    }) %>%
    lapply(., paste0, collapse = ";") %>%
    unlist()
}


for (col in c(
  "is_novel",
  "gene_name",
  "transcript_name",
  "class_code",
  "comparison",
  "method",
  "as_type"
)) {
  if (col == "as_type") {
    df[[col]] <- simplify(df[[col]], remove = c("JS", "JE"))
    next
  }

  df[[col]] <- simplify(df[[col]])
}

# TODO html link

write.xlsx(
  df %>% select(coordinates, everything()),
  overwrite = T,
  firstCol = T,
  colNames = T,
  asTable = T,
  file = as.character(opt$output),
  sheetName = "baltica_table"
)