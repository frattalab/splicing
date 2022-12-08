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
  library(stringr)
  library(reshape2)
  library(readr)
})

#' Compute and filter hits based on the difference in the genomic start and end
#'
#' @param query first set of range
#' @param subject optional, second set of ranges
#' @param max_start max_end absolute max difference at the start and
#' end coordinates,
#' respectively
#' @return overlapping ranges given the constrain
#' @export
filter_hits_by_diff <- function(query,
                                subject = NULL,
                                max_start = 2,
                                max_end = 2) {
  stopifnot(is(query, "GRanges"))
  if (is.null(subject)) {
    hits <- findOverlaps(query)
    subject <- query
  } else {
    stopifnot(is(subject, "GRanges"))
    hits <- findOverlaps(query, subject)
  }
  query <- query[queryHits(hits)]
  subject <- subject[subjectHits(hits)]
  start_dif <- abs(start(query) - start(subject))
  end_dif <- abs(end(query) - end(subject))
  hits <- hits[start_dif <= max_start & end_dif <= max_end]
  hits
}

#' Creates exon by transcript list (ex_by_tx) and remove items
#' with a single exon
#'
#' @param gtf GRange file loaded with rtracklayer::import
#' @return  a list of exon by transcripts, excluding single exon trascripts
#' @export
filter_multi_exon <- function(gtf) {
  stopifnot(is(gtf, "GRanges"))

  ex <- subset(gtf, type == "exon")
  stopifnot(length(ex) > 0)

  # discard single exons transcripts
  multi_ex <- table(ex$transcript_id) > 1
  ex <- subset(ex, mcols(ex)$transcript_id %in% names(multi_ex[multi_ex]))
  ex_tx <- split(ex, ex$transcript_id)
  ex_tx
}

#' Compute a set of introns from GTF files
#'
#' @param gtf_path path to to the gtf file
#' @param read_gtf function to read the gtf_file
#' @return a GRange that obj with the introns named by their parent trancripts
#' @export
get_introns <- function(ex_tx) {
  stopifnot(is(ex_tx, "List"))

  introns <- psetdiff(unlist(range(ex_tx), use.names = FALSE), ex_tx)
  introns <- unlist(introns)
  introns
}

aggregate_annotation <- function(gr) {
  stopifnot(is(gr, "GRanges"))

  equal_hits <- findOverlaps(
    gr,
    type = "equal"
  )
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

#' Fetch exons pairs for an intron
#'
#' @param introns_by_transcript list of introns by transcript
#' @return a data.frame with acceptor and donor exon number for an intron
#' @export
get_exon_number <- function(ex_tx) {
  stopifnot(is(ex_tx, "List"))

  exon_n_by_transcript <- lapply(ex_tx, function(x) embed(x$exon_number, 2))
  exon_number <- tidyr::unnest(
    tibble::enframe(exon_n_by_transcript),
    cols = "value"
  )

  exon_number
}

#' Get the most frequent variable from a numeric
#' @param .x numeric vector
#' @return the most frequent numeric element from .x
most_frequent <- function(.x) {
  as.numeric(names(which.max(table(.x))))
}

#' Process a file with scores from orthogonal experiment
#' Used to integrated third-generation sequencing results into baltica
#'
#' File formats should be readable by rtracklayer (BED, GTF)
#' the score column is mandatory
#' the comparison column is facultative
#' @param .path file path
#' @return a GRange object with a score metadata column
#' @export
process_ort_result <- function(.path) {
  message("Processing orthogonal result")

  ort_result <- rtracklayer::import(.path)
  if (!"score" %in% colnames(mcols(ort_result))) {
    warn(
      str_glue("score column is missing from {.path}
      , so the file will not be used")
    )
    ort_result <- GRanges()
  } else {
    .score <- mcols(ort_result)$score
    if ("comparison" %in% colnames(mcols(ort_result))) {
      .comparison <- mcols(ort_result)$comparison
    }
    mcols(ort_result) <- NULL
    mcols(ort_result)$score <- .score
    if (exists(".comparison")) {
      mcols(ort_result)$comparison <- .comparison
    }
  }
  ort_result
}

#' Corrects coordinates from alt GRange to ref GRange
#' by the most common start and end difference
#' @param ref reference GenomicRange
#' @param alt alternative GenomicRange
#' @return alt but corrected start and end coordinates
#' @export
integrate_coordinates <- function(ref, alt) {
  hits <- findOverlaps(ref, alt)
  query <- ref[queryHits(hits)]
  subject <- alt[subjectHits(hits)]
  start_diff <- start(query) - start(subject)
  start(alt) <- start(alt) + most_frequent(start_diff)
  end_diff <- end(query) - end(subject)
  end(alt) <- end(alt) + most_frequent(end_diff)
  return(alt)
}

default_input <- str_glue(
  "{method}/{method}_junctions.csv",
  method = c("junctionseq", "leafcutter", "majiq", "rmats")
) %>% str_c(collapse = ",")

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Path to parsed DJU result",
    metavar = "character",
    default = default_input
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
  ),
  make_option(
    c("-n", "--orthogonal_result"),
    type = "character",
    help = "Optionall, path to a file readable by rtracklayer (BED, GFF) with
     a score column with results from a orthogonal method, such as
     paired ONT Nanopore-seq from the same samples.",
    metavar = "character",
  ),
  make_option(
    c("--unstranded"),
    type = "logical",
    help = "If library is unstranded, [default is %default]",
    default = FALSE
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
    output = snakemake@output[[1]],
    orthogonal_result = snakemake@params$orthogonal_result,
    unstranded = snakemake@params$unstranded
  )
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
}

files <- strsplit(opt$input, ",")[[1]]

if (!all(as.logical(lapply(files, file.exists)))) {
  stop("Input file not found.", call. = FALSE)
} else if (!file.exists(opt$annotation)) {
  stop("Annotation not found.", call. = FALSE)
}

majiq_idx <- grep("majiq", files)
leafcutter_idx <- grep("leafcutter", files)
junctionseq_idx <- grep("junctionseq", files)
rmats_idx <- grep("rmats", files)

message("Processing DJU method output")
df <- list(
  majiq = read.csv(files[[majiq_idx]]),
  leafcutter = read.csv(files[[leafcutter_idx]]),
  junctionseq = read.csv(files[[junctionseq_idx]]),
  rmats = read.csv(files[[rmats_idx]])
)

df$junctionseq["score"] <- 1 - df$junctionseq$padjust
df$leafcutter["score"] <- 1 - df$leafcutter$p.adjust
df$majiq["score"] <- 1 - df$majiq$probability_non_changing
df$rmats["score"] <- 1 - df$rmats$FDR
comparison <- unique(unlist(lapply(df, function(x) unique(x$comparison))))

message("Processing de novo annotation")
gtf <- rtracklayer::import.gff2(opt$annotation)
tx <- subset(gtf, type == "transcript")
ex_tx <- filter_multi_exon(gtf)
introns <- get_introns(ex_tx)

reference <- rtracklayer::import.gff(opt$reference)
ref_ex_tx <- filter_multi_exon(reference)
ref_introns <- get_introns(ref_ex_tx)
ref_introns <- ref_introns[width(ref_introns) > 2, ]

introns$is_novel <- !(introns %in% ref_introns)

if (!is.null(opt$orthogonal_result)) {
  ort_result <- process_ort_result(opt$orthogonal_result)

  mcols(ort_result)["method"] <- "orthogonal"
  # resolves missing comparison column if there is a single column
  if (length(comparison) == 1 & !("comparison" %in% colnames(ort_result))) {
    mcols(ort_result)["comparison"] <- comparison[[1]]
  }
  df$orthogonal <- ort_result
}


message("Preparing annotation")
# TODO prepate annotation in another script, this does not need be done here
# It's super slow

exon_number <- get_exon_number(ex_tx)
txid_to_gene <- setNames(tx$gene_name, nm = tx$transcript_id)
txid_to_tx_name <- setNames(tx$cmp_ref, nm = tx$transcript_id)
txid_to_classcode <- setNames(tx$class_code, nm = tx$transcript_id)

mcols(introns)$gene_name <- txid_to_gene[names(introns)]
mcols(introns)$transcript_name <- txid_to_tx_name[names(introns)]
mcols(introns)$class_code <- txid_to_classcode[names(introns)]
mcols(introns)$exon_number <- exon_number$value
introns$exon_number <- apply(introns$exon_number, 1, paste, collapse = "-")

introns <- aggregate_annotation(introns)

message("Proceding with integration")

if (isTRUE(opt$unstranded)) {
  message("Unstranded library, removing strand information")
  strand(introns) <- "*"
  strand(ref_introns) <- "*"
  introns$coordinates <- as.character(introns)
}

gr <- lapply(df, function(x) {
  .gr <- GRanges(x)
  .gr <- .gr[width(.gr) > 2, ]
  if (!"comparison" %in% colnames(mcols(.gr))) {
    mcols(.gr)["comparison"] <- "NA"
  }
  mcols(.gr) <- mcols(.gr)[c("score", "comparison", "method")]

  if (isTRUE(opt$unstranded)) {
    strand(.gr) <- "*"
  }

  integrate_coordinates(ref_introns, .gr)
})

gr <- unlist(as(gr, "GRangesList"))


mcols(gr)$method <- names(gr)
sj <- mcols(gr)
sj$coordinates <- as.character(gr)
sj <- as.data.frame(sj)
sj <- dcast(
  sj,
  coordinates ~ method + comparison,
  value.var = "score",
  fun.aggregate = max
)

is.na(sj) <- sapply(sj, is.infinite)

annotation <- as.data.frame(mcols(introns))

sj <- merge(sj, annotation, by = "coordinates", all.x = TRUE, all.y = FALSE)

for (col in c(
  "class_code", "gene_name", "transcript_name", "exon_number"
)) {
  sj[[col]] <- as.character(lapply(sj[[col]], paste0, collapse = ";"))
}

sj$is_novel[is.na(sj$is_novel)] <- TRUE
sj$group <- NULL
sj <- sj[
  order(
    rowSums(
      sj[
        grepl(
          x = colnames(sj), pattern = "vs"
        )
      ],
      na.rm = T
    ),
    decreasing = TRUE
  ),
]


readr::write_csv(sj, opt$output)