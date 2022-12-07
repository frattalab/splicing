#' Compute and filter hits based on the fraction of overlap
#'
#' @param query first set of range
#' @param subject second set of ranges
#' @param cutoff for keeping the ranges
#' @return overlapping ranges given the contrain
#' @export

filter_hits_by_fraction <- function(query, subject, cutoff = 0.99) {
  stopifnot(is(query, "GRanges"))
  stopifnot(is(subject, "GRanges"))
  hits <- findOverlaps(query, subject)
  query <- query[queryHits(hits)]
  subject <- subject[subjectHits(hits)]
  overlap <- width(pintersect(query, subject)) / pmin(width(query), width(subject))
  hits[overlap >= cutoff]
}


#' Compute and filter hits based on the difference in the genomic start and end
#'
#' @param query first set of range
#' @param subject second set of ranges
#' @param max_start max_end absolute max difference at the start and end coordinates, respectively
#' @return overlapping ranges given the contrain
#' @export

filter_hits_by_diff <- function(query, subject, max_start = 2, max_end = 2) {
  stopifnot(is(query, "GRanges"))
  stopifnot(is(subject, "GRanges"))
  hits <- findOverlaps(query, subject)
  query <- query[queryHits(hits)]
  subject <- subject[subjectHits(hits)]
  start_dif <- abs(start(query) - start(subject))
  end_dif <- abs(end(query) - end(subject))
  hits <- hits[start_dif <= max_start & end_dif <= max_end]
  hits
}

#' Creates exon by transcript list (ex_by_tx) and remove items with a single exon
#'
#' @param gtf GRange file loaded with rtracklayer::import
#' @return  a list of exon by transcripts, excluding single exon trascripts 
#' @export

filter_multi_exon <- function(gtf) {
  stopifnot(is(gtf, "GRanges"))

  ex <- subset(gtf, type == 'exon')
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


#' Fetch exons pairs for an intron
#'
#' @param introns_by_transcript list of introns by transcript
#' @return a data.frame with acceptor and donor exon number for an intron
#' @export
get_exon_number <- function(ex_tx) {
  stopifnot(is(ex_tx, "List"))

  exon_n_by_transcript <- lapply(ex_tx, function(x) embed(x$exon_number, 2))
  exon_number <- tidyr::unnest(tibble::enframe(exon_n_by_transcript ), cols = 'value')

  exon_number
}


#' Annotated a set of features with overllaping gene_name
#'
#' @param gtf annotation with gene_names 
#' @param df 
#' @return a data.frame with acceptor and donor exon number for an intron
#' @export

annotate_gene <- function(gtf, df) {
  stopifnot(is(gtf, "GRanges"))
  stopifnot(is(df, "data.frame"))
  gr <- GRanges(df)
  tx <- subset(gtf, type == 'transcript')
  hits <- findOverlaps(gr, tx)
  hits_by <- split(subjectHits(hits), queryHits(hits))
  hits_by <- lapply(hits_by, function(x)  paste0(unique(mcols(tx)[x, 'gene_name']), collapse = ';'))
  df[as.numeric(names(hits_by)), 'gene_name'] <- as.character(hits_by)
  df
}

#' Giving a GRange genomic Annotated a set of features with overllaping gene_name
#'
#' @param introns_by_transcript list of introns by transcript
#' @return a data.frame with acceptor and donor exon number for an intron
#' @export

aggregate_metadata <- function(gr){
    stopifnot(is(gr, "GRanges"))
    
    equal_hits <- findOverlaps(gr, type='equal')
    meta <- data.frame(mcols(gr)[to(equal_hits), ])
    meta$index <- from(equal_hits)
    agg_metadata <- aggregate(. ~ index, meta, paste, collapse=';')
    gr <- gr[unique(agg_metadata$index)]
    mcols(gr) <- subset(agg_metadata, select = -c(index))
    gr

}

#' Append a string to the basename of a filepath
#' @param filepath path to be modified
#' @param to_append string to be added to the basename
#' @return a character
#' @export
#'
append_to_basename <- function(filepath, to_append='_nomatch'){
  sub("\\.([^\\.]*)$", paste0(to_append, "\\.\\1"), filepath)
}

#' Append a string to the basename of a filepath
#' @param junction.gr GRange of SJ
#' @param to_append GR of exons
#' @return alternative splice type for every SJ and Exon
#' @export
#'
as.type <- function(junction.gr, exon.gr) {
  stopifnot(is(junction.gr, "GRanges"))
  stopifnot(is(exon.gr, "GRanges"))


  j.strand <- as.character(strand(junction.gr)[1])
  j.start <- as.integer(start(junction.gr)[1])
  j.end <- as.integer(end(junction.gr)[1])
  e.start <- as.integer(start(exon.gr)[1])
  e.end <- as.integer(end(exon.gr)[1])

  is.annotated <- dplyr::case_when(
    j.strand == '+' & j.end == e.start ~ "JE",
    j.strand == '+' & j.start == e.end ~ "JS",
    j.strand == '-' & j.end == e.start ~ "JE",
    j.strand == '-' & j.start == e.end ~ "JS",
    TRUE ~ "NA"
  )

  if (is.annotated != "NA") {
    return(is.annotated)
  }

  if (j.strand == '+'){
    a <- e.start - j.start
    b <- j.end - e.end
    c <- e.end - j.start
    d = j.end - e.start
  } else if (j.strand == '-'){
    a = j.end - e.end
    b = e.start - j.start
    c = j.end - e.start
    d = e.end - j.start
  } else {
    message('AS type definition needs strand information')
    return("NA")
  }

  type <- dplyr::case_when(
    all(c(a > 0, b > 0, c > 0, d > 0)) ~ "ES",
    all(c(a < 0, b > 0, c > 0, d > 0)) ~ "A5SS",
    all(c(a > 0, b < 0, c > 0, d > 0)) ~ "A3SS",
    all(c(a > 0, b > 0, c > 0, d > 0)) ~ "EXITRON",
    TRUE ~ "NA")
  return(type)
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
process_RMATS_ass <- function(df, start='flankingES', end='shortES', type, FDR=0.05){
  
  df <- df %>%
    dplyr::filter(FDR < !!FDR)  %>% 
    dplyr::select(chr, !!start, !!end, strand, comparison, FDR, IncLevelDifference) %>% 
    dplyr::rename(c(start=!!start, end=!!end)) %>% 
    mutate(start=pmin(start, end), end=pmax(start, end))
  
  df <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
  df <- df[width(df) > 1,]
  df <- unique(df)
  mcols(df)['type'] <- type

  suppressWarnings(try(seqlevelsStyle(df) <- 'Ensembl'))
  
  
  df  
}


#' Process for exon skipping and intron retention rMATs output files 
#' @param df dataframe from rMATs 
#' @param start name of the start column
#' @param end name of the end column
#' @param type flag for splice junction type 
#' @param FDR is the FDR cutoff
#' @return GenomicRange of the selected SJ 
#' @export
#'
process_RMATS <- function(df, start, end, type, FDR=0.05) {
  
  df <- df %>%   
    dplyr::filter(FDR < !!FDR) %>% 
    dplyr::select(chr, !!start, !!end, strand, comparison, FDR, IncLevelDifference) %>% 
    dplyr::rename(c(start=!!start, end=!!end))
  
  df <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = T)
  df <- unique(df)
  df$type <- type
  suppressWarnings(try(seqlevelsStyle(df) <- 'Ensembl'))
  
  df  
}

  