#read in the bed files for a given group
#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library("optparse")


####end of the helper functions####
option_list = list(
    make_option(c("-l", "--scallopbed"), type="character", default=NULL,
                help="assembled transcripts from the transcriptome_assembly workflow", metavar="character"),
    make_option(c("-t", "--stringtiebed"), type="character", default=NULL,
                help="one of the deltaPSI annotated CSV files produced by this workflow", metavar="character"),
    make_option(c("-o", "--outputname"), type="character",
                help="output folder name - folder path for the combined events", metavar="character",default = "merged_cryptic_events/")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

scallopbed = opt$scallopbed
stringtiebed = opt$stringtiebed
outputname = opt$outputname

dir.create(outputname, showWarnings = FALSE)
comparison_name = basename(scallopbed)


scallopcsv = fread(paste0(scallopbed, ".csv"))
scallopbed = fread(scallopbed)
stringtiecsv = fread(paste0(stringtiebed, ".csv"))
stringtiebed = fread(stringtiebed)

both_bed = scallopbed |> 
    rbind(stringtiebed) |> 
    unique()

both_csv = scallopcsv |> 
    rbind(stringtiecsv) |> 
    unique()

outputname_bed = file.path(outputname,comparison_name)
outputname_csv = file.path(outputname,gsub(".bed",".csv",comparison_name))

fwrite(both_bed,outputname_bed,col.names = FALSE,sep = "\t")
fwrite(both_csv,outputname_csv)

