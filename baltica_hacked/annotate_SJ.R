#!/usr/bin/env Rscript
# Title     : annotate_SJ from many DJU methods
# Objective : To annotate SJ called DS by many DJU methods with gene and
# transcript level information
# Created by: Thiago Britto-Borges (thiago.brittoborges@uni-heidelberg.de)
# Created on: 30.04.20
suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(rtracklayer)
  library(readr)
  library(yaml)
  library(igraph)
})

# from https://stackoverflow.com/a/15373917/1694714
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Path to parsed DJU result",
    metavar = "character",
    default = "junctionseq/junctionseq_junctions.csv,
    leafcutter/leafcutter_junctions.csv,
    majiq/majiq_junctions.csv,
    rmats/rmats_junctions.csv"
  ),
  make_option(
    c("-a", "--annotation"),
    type = "character",
    help = "Path to annotation (GTF/GFF format)",
    metavar = "character",
    default = "stringtie/merged/merged.combined.gtf"
  ),
  make_option(
    c("-r", "--reference"),
    type = "character",
    help = "Path to reference annotation (GTF/GFF format)",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    help = "Path to output",
    metavar = "character",
    default = "results/SJ_annotated.csv"
  )
)
# enabling both baltica or snakemake input
if (exists("snakemake")) {
  opt <- list(
    input = paste0(
      c(
        snakemake@input[[1]],
        snakemake@input[[2]],
        snakemake@input[[3]],
        snakemake@input[[4]]
      ),
      collapse = ","
    ),
    reference = snakemake@params$ref,
    annotation = snakemake@input[[5]],
    output = snakemake@output[[1]]
  )
  snakemake@source("utils.R")
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
  source(file.path(dirname(thisFile()), "utils.R"))
}

files <- strsplit(opt$input, ",")[[1]]

message("Loading input")

if (!all(as.logical(lapply(files, file.exists)))) {
  stop("Input file not found.", call. = FALSE)
} else if (!file.exists(opt$annotation)) {
  stop("Annotation not found.", call. = FALSE)
}

gtf <- rtracklayer::import.gff2(opt$annotation)
.read_csv <- function(x) {
  suppressWarnings(
    readr::read_csv(x, col_types = readr::cols(
      chr = readr::col_character(),
      seqnames = readr::col_character()
    ))
  )
}
majiq_idx <- grep("majiq", files)
leafcutter_idx <- grep("leafcutter", files)
junctionseq_idx <- grep("junctionseq", files)
rmats_idx <- grep("rmats", files)

message("Processing GRanges")
df <- list(
  majiq = .read_csv(files[[majiq_idx]]),
  leafcutter = .read_csv(files[[leafcutter_idx]]),
  junctionseq = .read_csv(files[[junctionseq_idx]]),
  rmats = .read_csv(files[[rmats_idx]])
)

gr <- suppressWarnings(c(
  GRanges(df$majiq),
  GRanges(df$leafcutter),
  GRanges(df$junctionseq),
  GRanges(df$rmats)
))

mcols(gr) <- NULL
mcols(gr)$method <- c(
  df$majiq$method,
  df$leafcutter$method,
  df$junctionseq$method,
  df$rmats$method
)
mcols(gr)$comparison <- c(
  df$majiq$comparison,
  df$leafcutter$comparison,
  df$junctionseq$comparison,
  df$rmats$comparison
)
mcols(gr)$score <- c(
  df$majiq$probability_non_changing,
  df$leafcutter$p.adjust,
  df$junctionseq$padjust,
  df$rmats$FDR
)

message("Processing annotation")
tx <- subset(gtf, type == "transcript")
tx <- subsetByOverlaps(tx, gr)
gtf <- gtf[gtf$transcript_id %in% unique(tx$transcript_id)]

ex_tx <- filter_multi_exon(gtf)
introns <- get_introns(ex_tx)

if (!is.null(opt$reference)) {
  message("Processing reference annotation")

  reference <- rtracklayer::import.gff(opt$reference)
  ref_ex_tx <- filter_multi_exon(reference)

  ref_introns <- get_introns(ref_ex_tx)
  introns$is_novel <- !(introns %in% ref_introns)
}

message("Annotating set of introns")
exon_number <- get_exon_number(ex_tx)
txid_to_gene <- setNames(tx$gene_name, nm = tx$transcript_id)
txid_to_tx_name <- setNames(tx$cmp_ref, nm = tx$transcript_id)
txid_to_classcode <- setNames(tx$class_code, nm = tx$transcript_id)

mcols(introns)$gene_name <- txid_to_gene[names(introns)]
mcols(introns)$transcript_name <- txid_to_tx_name[names(introns)]
mcols(introns)$class_code <- txid_to_classcode[names(introns)]
mcols(introns)$exon_number <- exon_number$value
introns$exon_number <- apply(introns$exon_number, 1, paste, collapse = "-")
names(introns) <- NULL

aggregate_annotation <- function(gr) {
  stopifnot(is(gr, "GRanges"))

  equal_hits <- findOverlaps(gr, type = "equal")
  equal_hits <- as.data.frame(equal_hits)
  equal_hits <- igraph::graph_from_data_frame(equal_hits)
  equal_hits_groups <- stack(igraph::groups(igraph::clusters(equal_hits)))

  gr_ <- as.data.frame(mcols(gr))
  gr_$coordinates <- as.character(gr)
  gr_[as.numeric(equal_hits_groups$values), "group"] <- equal_hits_groups$ind
  gr_ <- aggregate(. ~ group, gr_, unique)
  gr <- GenomicRanges::GRanges(gr_$coordinates)
  mcols(gr) <- gr_
  gr
}

introns <- aggregate_annotation(introns)

hits <- filter_hits_by_diff(gr, introns, max_start = 2, max_end = 2)

x <- mcols(gr)
y <- mcols(introns)
x[from(hits), "hits"] <- to(hits)
x[is.na(x$hits), "hits"] <- -1
y["hits"] <- seq_len(nrow(y))

# all(mcols(gr)$score == xy$score)

message("Matching SJ to introns")
xy <- plyr::join(data.frame(x), data.frame(y))
mcols(gr) <- xy

gr <- as.data.frame(gr)

for (col in c("gene_name", "transcript_name", "class_code", "exon_number")) {
  gr[[col]] <- unlist(lapply(
    gr[[col]], function(x) {
      paste0(
        unique(x),
        collapse =  ";"
      )
    }
  ))
}


readr::write_csv(gr, opt$output)