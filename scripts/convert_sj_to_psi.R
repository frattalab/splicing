#!/usr/bin/env Rscript

library("dasper")
library(GenomicRanges)
library(tidyverse)
library(data.table)
library("optparse")


output_the_psi_files = function(sample_name,
  sample_file,
  gtf,
  output_folder,
  mincount = 1){

  output_filepath_raw = glue::glue("{output_folder}/{sample_name}_annotated.csv")
  output_filepath_normed = glue::glue("{output_folder}/{sample_name}_normalized_annotated.csv")

  ref =  GenomicFeatures::makeTxDbFromGFF(gtf,format = 'gtf')

  junctions <-
      junction_load(
          junction_paths = sample_file
      )

  junctions <- junction_annot(junctions, ref)

  junctions_filtered <- junction_filter(junctions,
      count_thresh = c(raw = mincount)) %>% junction_norm()

  annotated_clustered= as.data.table(junctions_filtered@rowRanges) %>%
    dplyr::select(seqnames,start,end,strand_junction,type,gene_id_junction,index,clusters)

  normed = as.data.table(SummarizedExperiment::assays(junctions_filtered)[["norm"]])
  normed$index = 1:nrow(normed)


  annotated_clustered_normed = annotated_clustered %>% left_join(normed,by = "index")
g
  setnames(annotated_clustered_normed,"count_1", sample_name)

  annotated_clustered_normed = annotated_clustered %>% left_join(normed,by = "index")

  setnames(annotated_clustered_normed,"count_1", sample_name)


  fwrite(annotated_clustered_normed,output_filepath_normed)

  print("File all written!")

}

option_list = list(
    make_option(c("-n", "--sample_name"), type="character", default=NULL,
                help="the final sample output name", metavar="character"),
    make_option(c("-f", "--sample_file"), type="character", default=NULL,
                help="the final sample file path", metavar="character"),
    make_option(c("-g", "--gtf"), type="character", default=NULL,
                help="GTF to annotated against", metavar="character"),
    make_option(c("-o", "--output_folder"), type="character", default="out.txt",
                help="output file name", metavar="character")
    make_option(c("-m", "--mincount"), type="character", default=1
                help="mincount to consider a junction", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


output_the_psi_files(opt$sample_name,opt$sample_file,opt$gtf,opt$output_folder)
