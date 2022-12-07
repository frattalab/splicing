#!/usr/bin/env Rscript
# Title : JunctionSeq analysis
# Objective : Call functions for JunctionSeq analysis in context of Baltica
# Created by: tbrittoborges
# Created on: 25/04/2018
suppressPackageStartupMessages(library(JunctionSeq))
options(stringsAsFactors = FALSE);

args<-commandArgs(TRUE)

decoder <- read.table(
  args[1], header = T, stringsAsFactors = F);
threads <- as.integer(args[3])

sample.files <- paste0(
  "junctionseq/mergedOutput/",
  decoder$sample.ID,
  "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
message("Starting runJunctionSeqAnalyses (", date(), ")");

jscs <- runJunctionSeqAnalyses(
  sample.files = sample.files,
  sample.names = decoder$sample.ID,
  condition = decoder$group.ID,
  nCores = threads,
  verbose = TRUE,
  debug.mode = TRUE,
  flat.gff.file = 'junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz',
  use.multigene.aggregates = TRUE,
  analysis.type = "junctionsOnly");

prefix <- sub('sigGenes.results.txt.gz', '', args[2])
message("Starting writeCompleteResults (", date(), ")")
writeCompleteResults(
  jscs,
  outfile.prefix = prefix,
  save.jscs = TRUE)

message("Done with writeCompleteResults (", date(), ")")
