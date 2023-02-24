#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(GenomicRanges)
  library(trackViewer)
  library(rtracklayer)
  library(UpSetR)
})

files <- Sys.glob('../stringtie/*/merged.gtf')
gtfs <- lapply(files, rtracklayer::import)
names(gtfs) <- str_split(files, '/', simplify = T)[, 3]

gtfs <- lapply(gtfs, subset, type == 'exon')

for (n in names(gtfs)){
  mcols(gtfs[[n]])['junction'] <- as.character(gtfs[[n]])
  gtfs[[n]] <- as_tibble(mcols(gtfs[[n]]))
}

gtfs <- bind_rows(gtfs, .id='id')
gtfs$cov <- as.numeric(gtfs$cov)

e <- gtfs %>%
  select(-c(source, type, score, phase, FPKM, TPM, gene_id, transcript_id, exon_number)) %>%
  group_by(junction, id) %>%
  summarise(cov_sum = sum(cov)) %>%
  ungroup() %>%
  group_by(junction) %>%
  spread(key = id, value = cov_sum) %>%
  ungroup()

e <- e[ rowSums(e[2: 9], na.rm = T) > 10, ]
e.na <- is.na(e[2:9])
e.na <- ifelse(e.na, 0, 1)

pdf('exons_upset.pdf')
upset(
  as.data.frame(e.na),
  nsets = 120,
  nintersects = 20,
  order.by = 'freq',
  scale.intersections = "log10", scale.sets = "log10",
  sets.x.label = "Number of Exons (log10)")
dev.off()

# work with the annotation

ref <- import('/biodb/genomes/mus_musculus/GRCm38_90/GRCm38.90.gtf')
ref <- subset(ref, type == 'exon')
gr <- GRanges(e)

olaps <- findOverlaps(ref, gr, type='equal', ignore.strand=TRUE)

e[subjectHits(olaps), 'annotated'] <- TRUE
e$annotated[ is.na(e$annotated) ] <- FALSE

e.na <- as.data.frame(e.na)
e.na$annotated <- e$annotated


is.annotated <- function(x) {
  x["annotated"] == TRUE
}

pdf('exons_upset_withannotation.pdf')
upset(
  as.data.frame(e.na),
  queries = list(list(query = is.annotated, params = list(), active = T)),
  nsets = 120,
  nintersects = 20,
  # scale.intersections = "log10", scale.sets = "log10",
  sets.x.label = "Number of Exons",
  mainbar.y.label = "Intersection Size",
  text.scale = 1.8,
  group.by = "degree",
  order.by="freq",
  mb.ratio = c(0.5, 0.5))
dev.off()




