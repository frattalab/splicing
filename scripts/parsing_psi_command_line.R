#!/usr/bin/env Rscript
library("optparse")
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="psi single tsv file", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt", 
                help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

parse_voila_tsv = function(file_path){
    library(data.table)
    library(splitstackshape)
    library(janitor)
    ####this function takes the file path to a voila tsv and parses it using splitstackshape and then combines the lsv_id and junction coordinates to make a unique
    ###lsv_id_junction for each lsv to get a PSI per splice event in long format
    #read in the file using data table and janitor to get clean names
    if(file.exists(file_path)){
        voila_org =  as.data.table(clean_names(fread(file_path)))
        colsToSplit = c("e_psi_per_lsv_junction",
                        "var_e_psi_per_lsv_junction", 
                        "junctions_coords",
                        "exons_coords",
                        "de_novo_junctions")
        voila_melt = cSplit(voila_org, colsToSplit,';', "long")
        #the splitting and melting can produce non unique rows, so remove those
        voila_melt = unique(voila_melt)
        #also can result in some rows that have NA in their psi columns
        voila_melt = voila_melt[!is.na(e_psi_per_lsv_junction)]
        # split the exon_coords to have an exon start and exon end
        voila_melt[,c("exon_start","exon_end") := tstrsplit(exons_coords, "-")]
        # same with the junctions and the retained introns
        #now to make a unique identifier cmbined the lsv_id and junction coordinates
        voila_melt[,lsv_junc := paste(lsv_id, junctions_coords,sep = "_")]
        
        voila_melt[,c("junc_start","junc_end") := tstrsplit(junctions_coords, "-")]
        voila_melt[,c("ir_start","ir_end") := tstrsplit(ir_coords, "-")]
        # I find this column irritating
        voila_melt$ucsc_lsv_link = NULL
        # remove these as they're now redundant
        voila_melt$junctions_coords = NULL
        voila_melt$ir_coords = NULL
        voila_melt$exons_coords = NULL
        
        # helpful column
        voila_melt[,paste_into_igv_junction := paste0(chr, ":",junc_start, "-",junc_end)]
        # add on the target or source
        voila_melt[,exon_type := tstrsplit(lsv_type,"|")[1]]
        voila_melt[,exon_type := ifelse(exon_type == "t",'target',"source")]
        
        # this is a useless column for us
        voila_melt$lsv_type = NULL
        return(voila_melt)
    } else{
        print(file_path)
        return("File not found")
    }
}



parsed = parse_voila_tsv(opt$input)
fwrite(parsed,opt$out)