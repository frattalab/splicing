#!/usr/bin/env Rscript
# Title     : annotate_SJ from many DJU methods
# Objective : To annotate SJ called DS by many DJU methods with gene and transcript level information
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
    help = "Path to parsed DJU result",,
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

opt <- parse_args(OptionParser(option_list = option_list))
df <- suppressMessages(read_csv(opt$input))

simplify <- function(x, remove=c()){   
    x  %>% str_split(';') %>%
    lapply(., function(x) { setdiff(x, remove) }) %>% 
    lapply(., paste0, collapse = ';')   %>% 
    unlist()
}


df <- df %>% 
    select(-c(width)) %>% 
    mutate(intron = as.character(str_glue('{seqnames}:{start}-{end}:{strand}'))) %>% 
    select(-c(seqnames, start, end, strand))


for (col in c('is_novel', 'gene_name', 'transcript_name', 'class_code', 'comparison', 'method', 'as_type')){
    if (col == 'as_type'){
        df[[col]]  <- simplify(df[[col]], remove = c('JS', 'JE'))    
        next
    }
    
    df[[col]]  <- simplify(df[[col]])
    
}

# TODO provide more expressive output for score and AS_Type
# for example: {method}_{comparison}_{score}
# {transcript_id}_{exon_id}_{AS_assign}
openxlsx::write.xlsx(
  df,
  file = opt$output,
  sheetName = "Sheet1",
)
