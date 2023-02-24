#! /usr/bin/env Rscript
## ---------------------------
## Script name: create_config.R
## Purpose of script: interatively create baltica configuration
## Author: Thiago Britto-Borges
## Date Created: 2021-05-17
## Copyright (c) Thiago Britto-Borges, DieterichLab 2021
## Email: thiago.brittoborges@uni-heidelberg.de
suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(yaml)
})

workdir='~/sirv_test/'
gtf='~/Baltica/data/SIRV.gtf'
stopifnot(file.exists(gtf))
fasta='~/Baltica/data/SIRV.fa'
stopifnot(file.exists(fasta))
setwd(workdir)

files=Sys.glob('~/Baltica/data/*.bam')
length(files)

base_path=file.path(getwd(), 'analysis', 'baltica')
dir.create(base_path, recursive = T)

groups = case_when(
    grepl('mix1', files) ~ "mix1",
    grepl('mix2', files) ~ "mix2",
    grepl('mix3', files) ~ "mix3")

ids=paste(
    groups,
    ave(rep(NA, length(groups)), groups, FUN=seq_along),
    sep = '_' 
)

groups

ids

contrasts=list(
    c('mix2', 'mix1'),
    c('mix3', 'mix1'),
    c('mix3', 'mix2'))
names(contrasts)=sapply(contrasts, paste0, collapse="-vs-")

sample_path=normalizePath(dirname(dirname(files[1])))
files=gsub(pattern = 'workflow/mapping', x=files, '')

length(files)

samples=as.list(setNames(nm=ids, files))

config=dplyr::lst(
    path=base_path,
    sample_path=sample_path,
    assembly='NONE',
    ref=gtf,
    ref_fa=fasta,
    strandness='reverse',
    read_len=100L,
    samples=samples,
    contrasts=contrasts,
    qc_env_prefix="source ~/miniconda3/etc/profile.d/conda.sh; conda init; conda activate qc;",
    leafcutter_min_samples_per_group=2L,
    leafcutter_min_samples_per_intron=2L
)                                       

config

yaml::write_yaml(config, file.path(base_path, 'config.yaml'))

yaml::write_yaml(
    list("__default__" = list(
            'mem' = 16000L, 
            'cpu' = 4L,
            'out' = 'logs/{rule}_{wildcards}.log',
            'partition'='general'),
        "rseqc_read_duplication" = list(
            'mem' = 64000L),
        "build" = list(
            'mem' = 6400L, 
            'cpu' = 4L,
            'partition'='long'),
        "qc" = list(
            'mem' = 32000L),
        "junctioseq_analysis" = list(
            'mem' = 64000L, 
            'cpu' = 10L),
        "run_rmats" = list(
            'mem' = 64000L)
    ),
    file.path(base_path, 'cluster.yaml')
)


