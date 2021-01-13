#!/usr/bin/env Rscript

library("optparse")
library(data.table)
library(tidyverse)
option_list = list(
    make_option(c("-f", "--folder"), type="character", default=NULL,
                help="folder with the parsed csvs", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt",
                help="output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


smoosh_psi = function(folder,output){
    suffix = "_normalized_annotated.csv"
    # everything has a suffix that we define, use this for removing later
    files = list.files(folder,full.names = TRUE,pattern = suffix)


    # using something I found on the Google's

    mydata = fread(files[1])
    for(f in 2:length(files)){
      t = fread(files[f])
      t = t %>% dplyr::select(-index,-clusters,-gene_id_junction)
      mydata = full_join(mydata, t, by = c("seqnames","start","end","strand_junction","type"))
      rm(t)
    }

    fwrite(psi,output_psi)
    fwrite(variance,output_var)

}


smoosh_psi(folder = opt$folder, output = opt$out)
