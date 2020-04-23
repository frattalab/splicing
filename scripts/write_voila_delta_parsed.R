#!/share/apps/R-3.6.1/bin/Rscript
# V0.001
## Usage write_voila_delta_parsed.R infile outfile
## Produces list of cells that pass QC filter

# args <- commandArgs(trailingOnly=TRUE)
source("parsing_voila_delta.R")

#!/usr/bin/env Rscript

library("optparse")
option_list = list(
    make_option(c("-d", "--deltafile"), type="character", default=NULL,
                help="a majiq delta psi file", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt",
                help="output file name", metavar="character"),
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


output = parse_voila_delta_tsv(opt$deltafile)

fwrite(output = opt$out)
