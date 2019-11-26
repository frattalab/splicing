

library("SGSeq")
# make a file path funciton wtihout that annoying double slash :)
file.path2 = function(..., fsep = .Platform$file.sep){
  gsub("//", "/", file.path(..., fsep = fsep))
}
gtf = snakemake@input[[1]]

output_dir = snakemake@params[[1]]
n_cores = snakemake@threads[[1]]
sgseq.anno = file.path2(output_dir, "sgseq.anno.RData")
print( "creating SGSeq transcript annotations")
print( paste( "from", gtf ) )


# tx <- importTranscripts(gtf)
# txf_ucsc <- convertToTxFeatures(tx)
# sgf_ucsc <- convertToSGFeatures(txf_ucsc)
# sgv <- findSGVariants(sgf_ucsc)

tx <- importTranscripts(gtf)
txf <- convertToTxFeatures(tx)
sgf <- convertToSGFeatures(txf)
sgv <- findSGVariants(sgf,cores = n_cores)

print( paste( "saving to:", sgseq.anno )  )
save(tx, txf, sgf, sgv, file = sgseq.anno)
