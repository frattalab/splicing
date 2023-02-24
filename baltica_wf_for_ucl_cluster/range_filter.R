#' Compute and filter hits based on the fraction of overlap
#'
#' @param query first set of range
#' @param subject second set of ranges
#' @param cutoff for keeping the ranges
#' @return overlapping ranges given the contrain
#' @export
filter_hits_by_fraction <- function(query, subject, cutoff=0.99){
  stopifnot(is(query, "GRanges"))
  stopifnot(is(subject, "GRanges"))
  hits <- findOverlaps(query, subject)
  query <- query[queryHits(hits)]
  subject <- subject[subjectHits(hits)]
  overlap <- width(pintersect(query , subject)) / pmin(width(query), width(subject))
  hits[overlap >= cutoff]
}

#' Compute and filter hits based on the difference in the genomic start and end
#'
#' @param query first set of range
#' @param subject second set of ranges
#' @param max_start, max_end absolute max difference at the start and end coordinates, respectively
#' @return overlapping ranges given the contrain
#' @export
filter_hits_by_diff <- function(query, subject, max_start=2, max_end=2){
  stopifnot(is(query, "GRanges"))
  stopifnot(is(subject, "GRanges"))
  hits <- findOverlaps(query, subject)
  query <- query[queryHits(hits)]
  subject <- subject[subjectHits(hits)]
  start_dif <- abs(start(query) - start(subject))
  end_dif <- abs(end(query) - end(subject))
  hits <- hits[start_dif <= max_start &  end_dif <= max_end ]
}

