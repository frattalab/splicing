#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyr)
  library(readr)
  library(yaml)
  library(stringr)
  library(AnnotationHub)
})

sec.last.item <- function(x) {
  rev(x)[2]
}

get.names <- function(x) {
  sec.last.item(str_split(x, '/')[[1]])
}
setwd('/prj/Andre_Schneider_BNauheim/raw_silke_kreher/newBatchMar2019/Baltica')

# JunctionSeq ####
js <-
  read_delim('junctionseq/analysis/sigGenes.results.txt.gz', delim = '\t')
dim(js)
# js <- filter(js, padjust < 0.05)
# dim(js)

ah <- AnnotationHub()

query(ah, c("OrgDb", "Mus musculus"))
mm <- ah[['AH66157']]
#js$symbol <-
# AnnotationDbi::mapIds(mm,
#                       as.character(js$geneID),
#                       column = "SYMBOL",
#                       keytype = "ENSEMBL")
js$j <- str_glue_data(js, "{chr}:{start}-{end}")

# Majiq ####
majiq.files <- Sys.glob('majiq/*/*.deltapsi.tsv')
majiq.df <- lapply(majiq.files, read_delim, delim = '\t')
names(majiq.df) <- lapply(str_split(majiq.files, '/'),
                          '[[', 2)
majiq.df <- bind_rows(majiq.df, .id = 'id')
colnames(majiq.df) <-  make.names(colnames(majiq.df))
majiq.df <- majiq.df %>%
  unite('coords', Junctions.coords, IR.coords, sep = ';')

majiq.df$coords <- str_remove(majiq.df$coords, ';NA')

majiq.df <- separate_rows(
  majiq.df,
  c(
    E.dPSI..per.LSV.junction,
    P..dPSI...0.20..per.LSV.junction,
    P..dPSI...0.05..per.LSV.junction,
    coords
  ),
  sep = ';'
)

majiq.df$P..dPSI...0.20..per.LSV.junction <-
  as.numeric(majiq.df$P..dPSI...0.20..per.LSV.junction)
majiq.df$E.dPSI..per.LSV.junction <-
  as.numeric(majiq.df$E.dPSI..per.LSV.junction)

majiq.df$start <- as.integer(lapply(str_split(majiq.df$coords, '-'), '[[', 1))
majiq.df$stop <- as.integer(lapply(str_split(majiq.df$coords, '-'), '[[', 2))

# majiq is missing the chr column
library(rtracklayer)
gff <- readGFF('majiq/ref.gff')
head(gff)

head(majiq.df$Gene.ID)
gff <- subset(gff, type == 'gene')
id.2.chr <- setNames(gff$seqid, gff$ID)
majiq.df$chr <- id.2.chr[ majiq.df$Gene.ID]
majiq.df$j <- str_glue_data(majiq.df, "{chr}:{start}-{stop}")


# LeafCutter ####
# parse cluster results
cs.files <- Sys.glob('leafcutter/*/*_cluster_significance.txt')
cs <- lapply(cs.files, read_tsv)
names(cs) <- lapply(cs.files, get.names)
cs <-  bind_rows(cs, .id = 'contrast')

# parse effect sizes results
effectsizes.files <- Sys.glob('leafcutter/*/*effect_sizes.txt')
es <- lapply(effectsizes.files, read_tsv)
names(es) <- lapply(effectsizes.files, get.names)
es <- bind_rows(es, .id = 'contrast')

es$intron <- str_replace_all(es$intron, '_', ':')
intron <- read_delim(
  es$intron,
  delim = ':',
  col_names = c('chr', 'start', 'end', 'clu', 'clu_number', 'strand')
)

es <- select(es,-intron)
es <- bind_cols(es, intron)
es$cluster <- str_glue_data(es, "{chr}:{clu}_{clu_number}_{strand}")
es <- select(es,-c(clu, clu_number))

lc.df <- inner_join(es, cs)
lc.df$chr <- gsub('chr', '', lc.df$chr)
lc.df$j <- str_glue_data(lc.df, "{chr}:{start}-{end}")


# Integration
library(UpSetR)

# Upset all junctions
UpSetR::upset(
  UpSetR::fromList(
    list(
      majiq = majiq.df$j,
      leafcutter = lc.df$j,
      junctionSeq = js$j
      )
    ), 
  empty.intersections = "on", order.by = "freq")

# Upset significant junctions
subset(js, padjust < 0.05, select=j)
UpSetR::upset(
  UpSetR::fromList(
    list(
      majiq = subset(majiq.df, P..dPSI...0.20..per.LSV.junction > .1, select=j),
      leafcutter = subset(lc.df, p.adjust < 0.05, select=j),
      junctionseq = subset(js, padjust < 0.05, select=j))
    ),
  empty.intersections = "on", order.by = "freq"
)
# no significant junctions 


dat <- data.frame(
  junction = common,
  leafcutter=subset(
    lc.df, (j %in% common & contrast == 'dKO-vs-KO'), select=c(p.adjust)),
  majiq=subset(
    majiq.df, (j %in% common & id == 'dKO-vs-KO'), select=c(P..dPSI...0.20..per.LSV.junction))
)


df <- inner_join(
  lc.df, majiq.df, by=c("j" = "j", "contrast" = "id"))

ggplot(
  df, aes(x=p, y=P..dPSI...0.20..per.LSV.junction, color=contrast)) +
  geom_point(alpha=.2) +
  xlab('LeafCutter pvalue') +
  ylab('Majiq post prob') + 
  geom_vline(aes(xintercept=0.05),
             color="red", linetype="dashed", size=.2) +
  geom_hline(aes(yintercept=0.95),
            color="red", linetype="dashed", size=.2)

bw <- 2 * IQR(df$p) / length(df$p)^(1/3)
ggplot(df, aes(x=p, fill=contrast)) +
  geom_density(alpha=.05)

ggplot(df, aes(x=p, fill=contrast)) +
  geom_histogram(binwidth=bw, alpha=.2, position="identity")



subset(js.df, p.adjust < .05)
table(subset(lc.df, p.adjust < .05, select=c(contrast)))
table(subset(majiq.df, P..dPSI...0.20..per.LSV.junction > .80, select=c(id)))

intersect(
  subset(majiq.df, P..dPSI...0.20..per.LSV.junction > .75, select=j),
  subset(lc.df, p < 0.05, select=j)
)

subset(
  df,
  (p < 0.05 & P..dPSI...0.20..per.LSV.junction > .75),
  select=j)
  
