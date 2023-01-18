#!/usr/bin/env Rscript

# Title     : TODO
# Objective : TODO
# Created by: thiago
# Created on: 10.12.19

library(rtracklayer)
gtf <- import('/biodb/genomes/mus_musculus/GRCm38_87/GRCm38.87.gtf')
gtf_coding <- subset(gtf, transcript_biotype ==  'protein_coding')
export(gtf_coding, 'GRCm38.87_coding.gtf')


