library(pacman)

p_load("data.table","janitor","splitstackshape","gdata","dplyr","tibble")

parse_voila_tsv = function(file_path){
  ####this function takes the file path to a voila tsv and parses it using splitstackshape and then combines the lsv_id and junction coordinates to make a unique
  ###lsv_id_junction for each lsv to get a PSI per splice event in long format
  #read in the file using data table and janitor to get clean names
  if(file.exists(file_path)){
    voila_org =  as.data.table(clean_names(fread(file_path)))
    voila_melt = cSplit(voila_org, c("e_psi_per_lsv_junction","var_e_psi_per_lsv_junction", "junctions_coords","de_novo_junctions"),';', "long")
    #the splitting and melting can produce non unique rows, so remove those
    voila_melt = unique(voila_melt)
    #also can result in some rows that have NA in their psi columns
    voila_melt = voila_melt[!is.na(e_psi_per_lsv_junction)]
    #now to make a unique identifier cmbined the lsv_id and junction coordinates
    voila_melt[,lsv_junc := paste(lsv_id, junctions_coords,sep = "_")]
    return(voila_melt)
  } else{
    print(file_path)
    return("File not found")
  }
}

named.list <- function(...) { 
  l <- list(...)
  names(l) <- as.character( match.call()[-1] )
  return(l)
}

nans_with_zero = function(DT) {
  
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
}

psi_voila = "/Users/annaleigh/Documents/data/ursa_work/psi_voila_tsv_single/"
suffix = ".Aligned.sorted.out.psi.tsv"



args <- commandArgs(trailingOnly=TRUE)
gtf <- fread(args[1])