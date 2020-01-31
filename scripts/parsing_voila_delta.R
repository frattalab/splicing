library(pacman)
p_load("data.table","janitor","splitstackshape")

parse_voila_delta_tsv = function(file_path){
  ####this function takes the file path to a voila tsv and parses it using splitstackshape and then combines the lsv_id and junction coordinates to make a unique
  ###lsv_id_junction for each lsv to get a PSI per splice event in long format
  #read in the file using data table and janitor to get clean names
  if(file.exists(file_path)){

    voila_org =  as.data.table(clean_names(fread(file_path)))
    # look for columns with conditional psi
    condition_psi_cols = colnames(voila_org)[grep("e_psi",colnames(voila_org))]
    # this changes depending on the flag we run to majiq
    changed_by_evidence = colnames(voila_org)[grep("p_d_psi_",colnames(voila_org))]
    columns_to_split = c("e_d_psi_per_lsv_junction",
                         changed_by_evidence,
                         "exons_coords",
                         "junctions_coords","de_novo_junctions",
                         "ir_coords",condition_psi_cols)
    voila_melt = cSplit(voila_org, splitCols = columns_to_split, sep = ';', direction = "long")
    #the splitting and melting can produce non unique rows, so remove those
    voila_melt = unique(voila_melt)
    #also can result in some rows that have NA in their psi columns
    voila_melt = voila_melt[!is.na(e_d_psi_per_lsv_junction)]
    #now to make a unique identifier cmbined the lsv_id and junction coordinates
    voila_melt[,lsv_junc := paste(lsv_id, junctions_coords,sep = "_")]
    # split the exon_coords to have an exon start and exon end
    voila_melt[,c("exon_start","exon_end") := tstrsplit(exons_coords, "-")]
    # same with the junctions and the retained introns
    voila_melt[,c("junc_start","junc_end") := tstrsplit(junctions_coords, "-")]
    voila_melt[,c("ir_start","ir_end") := tstrsplit(ir_coords, "-")]
    # I find this column irritating
    voila_melt$ucsc_lsv_link = NULL
    # remove these as they're now redundant
    voila_melt$junctions_coords = NULL
    voila_melt$ir_coords = NULL
    voila_melt$exons_coords = NULL
    # the 5% is the FDR
    setnames(voila_melt,"p_d_psi_0_05_per_lsv_junction", "FDR")
    return(voila_melt)
  } else{
    print(file_path)
    return("File not found")
  }
}

overlap_exons = function(voila_melt){
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
}

# file_path = "/Users/annaleigh/cluster/alb_projects/data/muscle/analysis/majiq/delta_psi_voila_tsv/notIbm_IBM.psi.tsv"
# voila_org =  as.data.table(clean_names(fread(file_path)))
# voila_org[,gene_simple := tstrsplit(gene_id, "\\.")[1]]
# 
# voila_melt = cSplit(voila_org, c("e_d_psi_per_lsv_junction","p_d_psi_0_10_per_lsv_junction", "exons_coords","junctions_coords","de_novo_junctions","not_ibm_e_psi", "ibm_e_psi"),';', "long")
# voila_melt = unique(voila_melt)
# voila_melt = voila_melt[!is.na(not_ibm_e_psi)]
# voila_melt[,lsv_junc := paste(lsv_id, junctions_coords,sep = "_")]
# voila_melt[,coords2 := paste0("chr",chr,":",junctions_coords)]
# 
# sig_voila_melt = voila_melt[p_d_psi_0_10_per_lsv_junction > 0.95]
# sig_voila_melt[,c("exon_start", "exon_end") := tstrsplit(exons_coords,"-")]
# 
# grange_sig = makeGRangesFromDataFrame(sig_voila_melt[!is.na(as.numeric(exon_start)) & !is.na(as.numeric(exon_end))], start.field = "exon_start", end.field = "exon_end", keep.extra.columns = TRUE,ignore.strand = FALSE)
# 
# temp = cbind(index = 1, exons = paste(unlist(subsetByOverlaps(exon_human,grange_sig[1])$exon_rank),collapse = ","))
# for(i in 2:length(grange_sig)){
#   temp = rbind(temp,cbind(index = i, exons = paste(unlist(subsetByOverlaps(exon_human,grange_sig[i])$exon_rank),collapse = ",")))
# }
