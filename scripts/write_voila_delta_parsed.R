#!/share/apps/R-3.6.1/bin/Rscript
# V0.001
## Usage write_voila_delta_parsed.R infile outfile
## Produces list of cells that pass QC filter

args <- commandArgs(trailingOnly=TRUE)
source("parsing_voila_delta.R")
