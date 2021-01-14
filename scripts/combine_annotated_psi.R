#!/usr/bin/env Rscript

library("optparse")
library(data.table)
library(tidyverse)
library(tidytext)
option_list = list(
    make_option(c("-f", "--folder"), type="character", default=NULL,
                help="folder with the parsed csvs", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt",
                help="output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# combine_annotated_psi = function(folder,output){
#     suffix = "_normalized_annotated.csv"
#     # everything has a suffix that we define, use this for removing later
#     files = list.files(folder,full.names = TRUE,pattern = suffix)
#
#
#     # using something I found on the Google's
#     pb = txtProgressBar(min = 2, max = length(files), initial = 2)
#     mydata = fread(files[1])
#     for(f in 2:length(files)){
#       t = fread(files[f])
#       t = t %>% dplyr::select(-index,-clusters,-gene_id_junction)
#       mydata = unique(full_join(mydata, t, by = c("seqnames","start","end","strand_junction","type")))
#       rm(t)
#       setTxtProgressBar(pb,f)
#     }
#
#     fwrite(mydata,output)
#
# }


# combine_annotated_psi(folder = opt$folder, output = opt$out)

combine_annotated_psi_long = function(folder,output){
    suffix = "_normalized_annotated.csv"
    # everything has a suffix that we define, use this for removing later
    files = list.files(folder,full.names = TRUE,pattern = suffix)


    # using something I found on the Google's

    mydata = fread(files[1])
    mydata$sample = colnames(mydata)[9]
    setnames(mydata,colnames(mydata)[9],"PSI")
    pb = txtProgressBar(min = 2, max = length(files), initial = 2)

    for(f in 2:length(files)){
      t = fread(files[f])
      t$sample = colnames(t)[9]
      setnames(t,colnames(t)[9],"PSI")
      mydata = rbind(mydata,t)
      rm(t)
            setTxtProgressBar(pb,f)

    }
    mydata[,paste_into_igv_junction := paste0(seqnames, ":",start, "-",end)]
    # my_sparse_data = cast_sparse(mydata, paste_into_igv_junction,sample, PSI)

    # saveRDS(my_sparse_data,output)
    fwrite(mydata,output)

}



combine_annotated_psi_long(folder = opt$folder, output = opt$out)
