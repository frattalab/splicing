#!/share/apps/R-3.6.1/bin/Rscript
# V0.001
## Usage write_voila_delta_parsed.R infile outfile

library("optparse")
library(data.table)
library(splitstackshape)
library(janitor)
parse_voila_delta_tsv = function(file_path){
  ####this function takes the file path to a voila tsv and parses it using splitstackshape and then combines the lsv_id and junction coordinates to make a unique
  ###lsv_id_junction for each lsv to get a PSI per splice event in long format
  #read in the file using data table and janitor to get clean names
  if(file.exists(file_path)){

    voila_org =  as.data.table(clean_names(fread(file_path)))
    # look for columns with conditional psi
    condition_psi_cols = colnames(voila_org)[grep("_mean_psi",colnames(voila_org))]


    columns_to_split = c("mean_dpsi_per_lsv_junction",
                        "probability_changing",
                         "probability_non_changing",
                         "junctions_coords",
                         "de_novo_junctions",
                         "ir_coords",condition_psi_cols)

    # for some reason....some of this columns didn't show up with me running majiq on the old bam files, different GFF or genome or something is possible route cuase sofor now
    voila_melt = cSplit(voila_org, splitCols = columns_to_split, sep = ';', direction = "long")
    #the splitting and melting can produce non unique rows, so remove those
    voila_melt = unique(voila_melt)
    #also can result in some rows that have NA in their psi columns
    voila_melt = voila_melt[!is.na(probability_changing)]
    #now to make a unique identifier cmbined the lsv_id and junction coordinates
    voila_melt[,lsv_junc := paste(lsv_id, junctions_coords,sep = "_")]
    # split the coords to have a start and an exon end
    # same with the junctions and the retained introns
    voila_melt[,c("junc_start","junc_end") := tstrsplit(junctions_coords, "-")]
    voila_melt[,c("ir_start","ir_end") := tstrsplit(ir_coords, "-")]
    voila_melt[,paste_into_igv_junction := paste0(seqid, ":",junc_start, "-",junc_end)]
    # I find this column irritating
    voila_melt$ucsc_lsv_link = NULL
    # remove these as they're now redundant
    voila_melt$junctions_coords = NULL
    voila_melt$ir_coords = NULL
    voila_melt$exons_coords = NULL
    # helpful columns
    voila_melt[,junc_dist := abs(as.numeric(junc_start) - as.numeric(junc_end))]



    return(voila_melt)
  } else{
    print(file_path)
    return("File not found")
  }
}


library("optparse")
option_list = list(
    make_option(c("-d", "--deltafile"), type="character", default=NULL,
                help="a majiq delta psi file", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt",
                help="output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


output = parse_voila_delta_tsv(opt$deltafile)
print('all done!')

fwrite(output, opt$out)
