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

sgfc_ucsc <- analyzeFeatures(si_quick_pm, features = txf_stand,cores = 3)

new_pred = annotate(sgfc_ucsc,txf_stand)

save(sgfc_ucsc, file = "~/ibm_pm_ctr.Rdata")

sgvc_new <- analyzeVariants(new_pred)

sgv_new <- rowRanges(sgvc_new)
sgvc_new <- getSGVariantCounts(sgv_new, sample_info = si_quick_pm)
sgvc

varcounts <- counts(sgvc_new)
vid <- variantID(sgvc_new)
eid <- eventID(sgvc_new)


# remove no count rows
noreads <- which(rowSums(varcounts) == 0 | rowSums(is.na(varcounts)) > 0)
if(length(noreads) > 0) {
  varcounts <- varcounts[-noreads,]
  vid <- vid[-noreads]
  eid <- eid[-noreads]
}


support.loc$condition <- factor(support[, condition])

# subset the samples from the matrix that are to be tested
varcounts.loc <- varcounts[, !is.na(support[,condition]) ]
message("varcounts subset OK")
#loc.countFiles <- countFiles[ !is.na(support.loc$condition) ]
support.loc <-  support.loc[ !is.na(support.loc$condition), ]
message("support subset OK")
# create formula
#formuladispersion <- ~ sample + (condition + type) * exon
#formula0 <-  ~ sample + condition
#formula1 <-  ~ sample + condition * exon



if( add.covariate == TRUE){
  names(support.loc)[ which( names(support.loc) == list.covariates) ] <- "type"
  formulaFullModel = ~ sample + exon + type:exon + condition:exon
  formulaReducedModel = ~ sample + exon + type:exon
}else{
  formulaFullModel = ~ sample + exon + condition:exon
  formulaReducedModel = ~ sample + exon
}

formulaFullModel = ~ sample_name + exon + condition:exon
formulaReducedModel = ~ sample_name + exon

sampleData = data.frame(sample_name = si_quick_pm$sample_name, condition = ifelse(grepl("IBM_",si_quick_pm$sample_name), "IBM", "not_IBM"))

DexSeqExons.loc <- DEXSeqDataSet(countData = varcounts,
                                 sampleData = sampleData,
                                 design = formulaFullModel,
                                 featureID = as.factor(vid),
                                 groupID = as.factor(eid))

DexSeqExons.loc <- estimateSizeFactors(DexSeqExons.loc)
DexSeqExons.loc <- DEXSeq::estimateDispersions(DexSeqExons.loc, formula = formulaFullModel)
DexSeqExons.loc <- DEXSeq::testForDEU(DexSeqExons.loc,
                                      reducedModel = formulaReducedModel,
                                      fullModel = formulaFullModel)
DexSeqExons.loc <- DEXSeq::estimateExonFoldChanges(DexSeqExons.loc, fitExpToVar="condition")
save(DexSeqExons.loc, file = "~/DexSeqExons.R")
res <- DEXSeq::DEXSeqResults (DexSeqExons.loc)
res.clean <- as.data.frame(res)

sample.data <- colData(sgvc_new)
sample.names <- sample.data$sample_name
n.samples <- length(sample.names)

countData.cols <- which( grepl("countData", names(res.clean) ) )
names(res.clean)[countData.cols] <- sample.names

res.clean$genomicData <- NULL
res.clean$FDR <- p.adjust(res.clean$pvalue, method = 'fdr')

sgvc.df <- as.data.frame(mcols(sgvc_new) )
psi <- variantFreq(sgvc_new)

if (length(noreads) > 0) {
  sgvc.df <- sgvc.df[-noreads,]
  psi <- psi[-noreads,]
}
# select just the samples present in the condition
psi <- psi[, sample.names]

colnames(psi) <- paste0(colnames(psi)  , "_psi")
res.clean <- cbind(res.clean, psi)
res.clean$geneName <- sgvc.df$geneName
res.clean$txName <- sgvc.df$txName
res.clean$eventID <- sgvc.df$eventID
res.clean$variantType <- sgvc.df$variantType
res.clean$from <- sgvc.df$from
res.clean$to <- sgvc.df$to
res.clean$type <- sgvc.df$type

annotations = convertGeneListTypes(unlist(res.clean$geneName),from_type = "ENTREZID",to_type = "ENSEMBL","hsapiens")
res.clean$external_gene_ID <- annotations$ENSEMBL[ match( res.clean$geneName, annotations$ENTREZID)]
res.clean$symbol <- annotations$SYMBOL[ match( res.clean$geneName, annotations$ENTREZID)]

createCoords <- function(from, to){
  fromSplit <- str_split_fixed(from, ":", 4)
  toSplit <- str_split_fixed(to, ":", 4)
  coord <- ifelse( toSplit[,4] == "+",
                   yes = paste0( toSplit[,2], ":", fromSplit[,3], "-", toSplit[,3]),
                   no = paste0( toSplit[,2], ":", toSplit[,3], "-", fromSplit[,3])
  )
}

res.clean$coords <- createCoords( from = res.clean$from, to = res.clean$to)
