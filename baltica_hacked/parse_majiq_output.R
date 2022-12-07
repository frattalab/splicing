#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(readr)
  library(dplyr)
  library(optparse)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "majiq/voila/*_voila.tsv",
    help = "Path with glob character to Majiq result files. [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "majiq/majiq_junctions.csv",
    help = "Path to output file [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cutoff"),
    default = 0.05,
    type = "double",
    help = "Discard junction with probability_non_changing < than --cutoff [default %default]"
  )
)
# enabling both baltica or snakemake input
if (exists('snakemake')) {
    opt <- list(
      input = snakemake@input,
      output = snakemake@output[[1]],
      cutoff = snakemake@params[['cutoff']]
      )
    files <- opt$input
    file_names <- gsub(
     x = opt$input,
     pattern = 'majiq/voila/(.+)_voila.tsv',
     replacement = '\\1')

  } else {

    opt <- parse_args(OptionParser(option_list = option_list))
    files <- Sys.glob(opt$input)
    file_names <- gsub(
      x = files,
      replacement = '\\1',
      pattern = sub(x=opt$input, pattern='\\*', replacement = '(.*)')
      )
  }

# rename column names from majiq result due to the presence of spaces
read_majiq_out <- function(x) {
  #           "Gene_Name Gene_ID LSV_ID E_dPSI_per_LSV_junction P_dPSI_beq_per_LSV_junction P_dPSI_leq_per_LSV_junction ref_E_PSI alt_E_PSI LSV_Type A5SS A3SS ES Num_Junctions Num_Exons De_Novo_Junctions chr strand Junctions_coords Exons_coords IR_coords UCSC_LSV_Link"  
  col_names <- "gene_name gene_id lsv_id mean_dpsi_per_lsv_junction probability_changing probability_non_changing ref_mean_psi alt_mean_psi lsv_type num_junctions num_exons de_novo_junctions chr strand junctions_coords exons_coords ir_coords ucsc_lsv_link"
  tmp <- read_tsv(
    x,
    col_names = strsplit(col_names, ' ')[[1]],
    comment = '#',
    cols(
      .default = col_character()
    )
  )
  tmp[-1,] # the skip parameter is not been registed in read_tsv 
}

res <- lapply(files, read_majiq_out)

names(res) <- file_names

res <- bind_rows(res, .id = 'comparison') %>%
  separate_rows(
    "junctions_coords",
    "mean_dpsi_per_lsv_junction",
    "probability_changing",
    "probability_non_changing",
    "ref_mean_psi",
    "alt_mean_psi",
    sep = ';',
    convert = T
  )
# flag canonical SJ
# res <- res %>%
#   arrange(comparison, lsv_id, ref_mean_psi) %>%
#   group_by(comparison, lsv_id) %>%
#   mutate(is_canonical = row_number() == 1) %>%
#   ungroup()


junction_pattern <- "(\\d+)-(\\d+)"
junctions_coords <- str_match(
  res$junctions_coords, junction_pattern)[, c(2, 3)]

res['start'] <- junctions_coords[, 1]
res['end'] <- junctions_coords[, 2]

message('Number of junctions output by Majiq ', nrow(res))
res <- dplyr::select(res, c(
  chr,
  start,
  end,
  strand,
  comparison,
  probability_changing,
  probability_non_changing,
  lsv_id,
  ref_mean_psi,
  alt_mean_psi)
) %>%
  filter(probability_non_changing < opt$cutoff) %>%
  mutate(method = 'majiq')

message('Number of junctions after filtering ', nrow(res))

write_csv(res, opt$output)
