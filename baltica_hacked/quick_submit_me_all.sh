#!/bin/bash

#Run calculate_splice_stability.py over all bam files in given folder
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=8G,h_rt=06:00:00,tmem=8G

# join stdout and stderr output
#$ -j y



cd /SAN/vyplab/alb_projects/pipelines/baltica_hacked2

source ~/.bash_profile

snakemake -s analysis.smk -c2 --configfile /SAN/vyplab/first_weeks/TDP_CHX_CLONES_GLIA/curves/baltica/config.yml --use-singularity --singularity-args "-B /SAN/vyplab:/SAN/vyplab"
