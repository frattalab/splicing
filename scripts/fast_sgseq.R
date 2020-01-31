library(stringr)
library(data.table)
library(SGSeq)
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

sample_name = gsub(".Log.final.out","",list.files(bam_dir,pattern = ".Log.final.out"))
si = data.table(sample_name = sample_name)
si[,file_bam := paste0(bam_dir,sample_name,bam_suffix,".bam"), by = 1:nrow(si)]
si[,paired_end := TRUE]
si[,read_length := 150]
si[, lib_size := getLibSize(paste0(bam_dir,sample_name,".Log.final.out")), by = 1:nrow(si)]

si[,frag_length := getInsertSize(paste0(bam_dir,sample_name,bam_suffix,".bam",".insert_size_metrics.txt")), by = 1:nrow(si)]

renaming = fread("/home/annbrown/muscle/rename.csv")

# renaming = fread('/Users/annaleigh/Documents/GitHub/muscle/rename.csv')
setDT(si)[renaming, new := type, on = c(sample_name = "current")]
si_quick_pm = si[new %like% "IBM|PM_|Ct"]
setnames(si_quick_pm, c("sample_name","new"), c("drop","sample_name"))
si_quick_pm$drop = NULL

library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txf_ucsc <- convertToTxFeatures(txdb)
seqlevelsStyle(txf_ucsc) = "NCBI"

txf_stand = keepStandardChromosomes(txf_ucsc,pruning.mode = "coarse")

sgfc_ucsc <- analyzeFeatures(si_quick_pm, features = txf_stand,verbose = TRUE, cores = 4)

save(sgfc_ucsc, file = "/home/annbrown/muscle/ibm_pm_ctr.Rdata")
