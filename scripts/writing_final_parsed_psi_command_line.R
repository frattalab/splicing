#!/usr/bin/env Rscript

library("optparse")
option_list = list(
    make_option(c("-f", "--folder"), type="character", default=NULL,
                help="folder with the parsed csvs", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt",
                help="output file name", metavar="character"),
    make_option(c("-s", "--suffix"), type="character", default=".Aligned.sorted.out",
                help="the suffix that everything gets appended with [default= '.Aligned.sorted.out_parsed.']", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


parse_all_the_parsed = function(folder,suffix,output){
    library(data.table)
    library(tidyverse)
    # get all the files that end in out_parsed_csv
    test_files = list.files(folder,full.names = TRUE)

    # everything has a suffix that we define, use this for removing later
    suffix = paste0(suffix,  "_parsed.csv")
    files = list.files(folder,full.names = TRUE,pattern = "_parsed.csv")


    # using something I found on the Google's
    mydata = tibble(File = files) %>%
        mutate(Data = lapply(File, read_csv)) %>%
        unnest(Data) %>%
        mutate(sample = stringr::str_split_fixed(File, "psi_voila_tsv_single/", 2)[,2]) %>%
        mutate(sample = gsub(suffix,"",sample)) %>%
        dplyr::select(-File)
    # doing this as 2 tables because I want to keep variance, but I don't
    # quite know how I want to use it yet
    psi = mydata %>%
        dplyr::rename(PSI = e_psi_per_lsv_junction) %>%
        dplyr::select(lsv_junc,chr,junc_start,junc_end,strand,PSI,sample) %>%
        pivot_wider(values_from = c("PSI"),names_from = 'sample')

    variance = mydata %>%
        dplyr::rename(varPSI = var_e_psi_per_lsv_junction) %>%
        dplyr::select(lsv_junc,chr,junc_start,junc_end,strand,varPSI,sample) %>%
        pivot_wider(values_from = c("varPSI"),names_from = 'sample')

    output_psi = file.path(output,"full_PSI.csv")
    output_var = file.path(output,"full_varPSI.csv")

    fwrite(psi,output_psi)
    fwrite(variance,output_var)

}


parse_all_the_parsed(folder = opt$folder, suffix = opt$suffix, output = opt$out)
