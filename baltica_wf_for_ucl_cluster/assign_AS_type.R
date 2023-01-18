#!/usr/bin/env Rscript
# Title     : assign_AS_type.R
# Objective : Assign AS type for events in Baltica table
# Created by: Thiago Britto-Borges (thiago.brittoborges@uni-heidelberg.de)
# Created on: 30.04.20
suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(reshape2)
  library(readr)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Path to parsed DJU result",
    metavar = "character",
    default = "results/SJ_annotated.csv"
  ),
  make_option(
    c("-a", "--annotation"),
    type = "character",
    help = "Path to annotation, in the GTF/GFF format",
    metavar = "character",
    default = "stringtie/merged/merged.combined.gtf"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    help = "Path to the output file",
    metavar = "character",
    default = "results/SJ_annotated_assigned.csv"
  )
)


if (exists("snakemake")) {
  opt <- list(
    input = snakemake@input[[1]],
    annotation = snakemake@input[[2]],
    output = snakemake@output[[1]]
  )
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
}


as.type <- function(junction.gr, exon.gr) {
  if (any(as.logical(strand(junction.gr) != strand(exon.gr)))) {
    warning("Junction and exon have different strand")
    return("NA")
  }

  if (length(junction.gr) != 1 | length(exon.gr) != 1) {
    warning("as.type function take one junction and exon per function call")
    return("NA")
  }

  j.strand <- as.character(strand(junction.gr)[1])
  j.start <- as.integer(start(junction.gr)[1])
  j.end <- as.integer(end(junction.gr)[1])
  e.start <- as.integer(start(exon.gr)[1])
  e.end <- as.integer(end(exon.gr)[1])

  is.annotated <- dplyr::case_when(
    j.strand == "+" & j.end == e.start ~ "JE",
    j.strand == "+" & j.start == e.end ~ "JS",
    j.strand == "-" & j.end == e.start ~ "JE",
    j.strand == "-" & j.start == e.end ~ "JS",
    TRUE ~ "NA"
  )

  if (is.annotated != "NA") {
    return(is.annotated)
  }

  if (j.strand == "+") {
    a <- e.start - j.start
    b <- j.end - e.end
    c <- e.end - j.start
    d <- j.end - e.start
  } else if (j.strand == "-") {
    a <- j.end - e.end
    b <- e.start - j.start
    c <- j.end - e.start
    d <- e.end - j.start
  } else {
    message("AS type definition needs strand information")
    return("NA")
  }

  type <- dplyr::case_when(
    all(c(a > 0, b > 0, c > 0, d > 0)) ~ "ES",
    all(c(a > 0, b < 0, c > 0, d > 0)) ~ "A5SS",
    all(c(a < 0, b < 0, c > 0, d > 0)) ~ "A3SS",
    TRUE ~ "NA"
  )
  return(type)
}

message("Loading input")

if (!file.exists(opt$input)) {
  stop("Input file not found.", call. = FALSE)
} else if (!file.exists(opt$annotation)) {
  stop("Annotation not found.", call. = FALSE)
}

gtf <- rtracklayer::import.gff2(opt$annotation)
df <- read_csv(opt$input)

message("Assigning AS type")
gr <- GRanges(df$coordinates)
exons <- subset(gtf, type == "exon")
hits <- as.data.frame(findOverlaps(gr, exons))
hits$as_type <- apply(hits, 1, function(x) as.type(gr[x[1]], exons[x[2]]))
.as_type <- stack(
  lapply(
    split(
      hits$as_type, hits$queryHits
    ), paste0,
    collapse = ";"
  )
)
.as_type <- rename(.as_type, values = "as_type")
df <- merge(x = df, y = .as_type, by.x = 0, by.y = "ind", all.x = T, sort = F)

message("Writting output")
df <- subset(df, select = -Row.names)
readr::write_csv(df, opt$output)