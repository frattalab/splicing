#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(optparse)
  library(Rsamtools)
})

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Differently spliced intron file", metavar="character"),
  # make_option(c("-o", "--out"), type="character", default=NULL,
  #             help="output file name [default= %default]", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL,
              help="Annotation in the GTF/GFF format", metavar="character"),
  make_option(c("-s", "--sequence"), type="character", default=NULL,
              help="Genomic sequence in the FASTA format", metavar="character")

);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (!file.exists(opt$input)){
  print_help(opt_parser)
  stop("Pass a valid intron file to -i or --input", call.=FALSE)
}
if (!file.exists(opt$annotation)) {
  print_help(opt_parser)
  stop("Pass a valid annotation to -a or --annotation ", call.=FALSE)
}
if (!file.exists(opt$sequence)) {
  print_help(opt_parser)
  stop("Pass a valid annotation to -a or --annotation ", call.=FALSE)
}


select.exons <- function(pos, exons, flank.dist=150) {
  dist <- distanceToNearest(pos, exons)

  if (any(is.na(subjectHits(dist)))) {
    warning('While trying to find nearest feature, missing values were created')
  }

  e <- exons[subjectHits(dist), ]
  genes <- unique(mcols(e)$gene_id)
  annotated <- mcols(dist)$distance == 0
  # for novel exons, use junction start + flank.dist nt for as the novel exon
  novel <- pos[queryHits(dist), ][!annotated, ]
  novel <- flank(novel, flank.dist, start = FALSE)

  mcols(novel) <- mcols(e[!annotated])
  e[!annotated] <- novel

  names(e) <- mcols(e)$gene_id
  score(e) <- 0

  exons.in.genes <- exons[mcols(exons)$gene_id %in% genes]
  e.negative <- GenomicRanges::setdiff(exons.in.genes, e)

  list(
    positive=e,
    negative=e.negative
  )
}

message('Loading input files.')
intron <- GRanges(read.csv(opt$input)) # 'Baltica/leafcutter/leafcutter_junctions_annotated.csv'
gtf <- rtracklayer::import.gff(opt$annotation) #  '/biodb/genomes/mus_musculus/GRCm38_90/GRCm38.90.gtf'
seq <- FaFile(opt$sequence)
intron <- split(intron, intron$contrast)
names(intron) <- gsub(pattern = '\\.', replacement='_', x= make.names(names(intron)))

exons <- gtf[gtf$type == 'exon']
exons <- exons[!duplicated(exons)]
# sequence.length <- quantile(width(exons), 0.5)

for (n in names(intron)) {
  message("Computing exons for the 5' splice site for ", n)
  x <- intron[[n]]

  fivepss <- ifelse(strand(x) == '+', end(x), start(x))

  acceptor <- GRanges(seqnames = seqnames(x),
                      IRanges(fivepss, width=1),
                      strand = strand(x))

  acceptor.exons <- select.exons(acceptor, exons)
  # rtracklayer::export(acceptor.exons[['positive']],
  #                     paste(n, 'acceptor_positive_exons.bed', sep = '_'))
  # rtracklayer::export(acceptor.exons[['negative']],
  #                     paste(n, 'acceptor_negative_exons.bed', sep='_'))

  writeXStringSet(
    getSeq(x = seq, acceptor.exons[['positive']]),
    paste(n, 'acceptor_positive_exons.fa', sep='_'))

  writeXStringSet(
    getSeq(x = seq, acceptor.exons[['negative']]),
    paste(n, 'acceptor_negative_exons.fa', sep='_'))


  message("Computing exons for the 3' splice site ", n)
  threepss <- ifelse(strand(x) == '+', start(x), end(x))
  donor <- GRanges(seqnames = seqnames(x),
                   IRanges(threepss, width=1),
                   strand = strand(x))

  donor.exons <- select.exons(donor, exons)
  # rtracklayer::export(donor.exons[['positive']],
  #                     paste(n, 'donor_positive_exons.bed', sep="_"))
  # rtracklayer::export(donor.exons[['negative']],
  #                     paste(n, 'donor_negative_exons.bed', sep="_"))

  writeXStringSet(
    getSeq(x = seq, donor.exons[['positive']]),
    paste(n, 'donor_positive_exons.fa', sep='_'))

  writeXStringSet(
    getSeq(x = seq, donor.exons[['negative']]),
    paste(n, 'donor_negative_exons.fa', sep='_'))



}
#
# for f in `ls *positive_exons.fa`
# do
# n=${f//positive/negative}
# sbatch --mem 32GB --partition long dreme-py3 -p $f -n $n -dna -norc -oc dreme/$f
# done
#
