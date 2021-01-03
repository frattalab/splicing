make_bed_from_annotated <- function(parsed_file, 
                                    output_filepath,
                                    only_up_in_case = FALSE,
                                    cutoff = NULL,
                                    junction_types = c("annotated", "none", "novel_acceptor", "novel_exon_skip", "novel_donor", 
                                                       "ambig_gene", "novel_combo")
                                    trackname = "All junctions"){
    #make the bed header
    bed_header = glue::glue('track name="{trackname}" description="{trackname}" visibility=2 itemRgb="On"')
    #read in the file, 
    annotated_junctions = data.table::fread(parsed_file)
    #keep only the junction types in the junction_types arguemnt
    if(only_up_in_case){
        #keep only the junctions increasing the case condition
        annotated_junctions = annotated_junctions[deltaPSI >0]
    }
    annotated_junctions = annotated_junctions[junc_cat %in% junction_types]
    #keep only the junctions with a abs(deltaPSI greater than cutoff)
    if(!is.null(cutoff)){
        annotated_junctions = annotated_junctions[abs(deltaPSI) > cutoff]
    }
    #make a column of gene name and the junction type
    annotated_junctions[,name := paste0(gene_name,":",junc_cat)]
    #keep only the relevant things
    annotated_junctions = unique(annotated_junctions[,.(seqnames,start,end,name,deltaPSI,strand,start,end)])
    #add colors
    annotated_junctions[,color := ifelse(deltaPSI > 0, "255,0,0","0,0,255")]
    writeLines(bed_header, output_filepath)
    write.table(annotated_junctions, output_filepath, 
                col.names=FALSE, 
                row.names=FALSE, 
                append=TRUE,sep = "\t",
                quote = FALSE)
    
    
)
    
}

####end of the helper functions####
option_list = list(
    make_option(c("-p", "--parsed"), type="character", default=NULL,
                help="parsed deltaPSI table produced by parse_voila_delta_tsv script", metavar="character"),
    make_option(c("-o", "--out"), type="character",
                help="output file name - no extension 2 files will be written", metavar="character",default = "annotated_junctions_bedformat.junctions.bed"),
    make_option(c("-c", "--cutoff"), type="character", default=NULL,
                help="deltaPSI cutoff for writing", metavar="character")
    make_option(c("-t", "--trackname"), type="character", default="All junctions",
                help="name of the track on IGV", metavar="character")
    make_option(c("-u", "--upincase"), type="character", default="All junctions",
                help="name of the track on IGV", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

make_bed_from_annotated(parsed_file = opt$parsed, 
                        output_filepath = opt$out, 
                        cutoff = opt$cutoff,
                        trackname = opt$trackname,
                        only_up_in_case = opt$upincase)