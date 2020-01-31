getLibSize <- function(logFile){

  stopifnot(file.exists(logFile))
  log <- readLines(logFile)
  unique <- log[ grepl("Uniquely mapped reads number", log)]
  multi <- log[ grepl("Number of reads mapped to multiple loci", log)]

  num_unique <- str_trim( str_split_fixed( unique, "\\|\t", 2)[,2] )
  num_multi <- str_trim( str_split_fixed( multi, "\\|\t", 2)[,2] )
  libSize <- as.numeric(num_unique) + as.numeric(num_multi)
  return(libSize)
}

getInsertSize <- function(logFile){

  if(file.exists(logFile)){
    log <- fread(logFile)$MODE_INSERT_SIZE

  }else(
    log = NA
  )


  return(log)
}
bam_dir = "/SAN/vyplab/alb_projects/data/muscle/analysis/STAR_aligned/"

# bam_dir = "/Users/annaleigh/cluster/alb_projects/data/muscle/analysis/STAR_aligned/"
bam_suffix = ".Aligned.sorted.out"
# make the sample information table
sample_name = gsub(".Log.final.out","",list.files(bam_dir,pattern = ".Log.final.out"))
si = data.table(sample_name = sample_name)
si[,file_bam := paste0(bam_dir,sample_name,bam_suffix,".bam"), by = 1:nrow(si)]
si[,paired_end := TRUE]
si[,read_length := 150]
si[, lib_size := getLibSize(paste0(bam_dir,sample_name,".Log.final.out")), by = 1:nrow(si)]

si[,frag_length := getInsertSize(paste0(bam_dir,sample_name,bam_suffix,".bam",".insert_size_metrics.txt")), by = 1:nrow(si)]
