#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y


SINGULARITY_CACHEDIR=/SAN/vyplab/alb_projects/tools/
export SINGULARITY_CACHEDIR


snakemake -s junctionseq.smk \
--configfile /SAN/vyplab/first_weeks/TDP_CHX_CLONES_GLIA/curves/baltica/config.yml \
--cluster-config /SAN/vyplab/alb_projects/pipelines/baltica_hacked2/cluster.yaml \
--jobscript /SAN/vyplab/alb_projects/pipelines/baltica_hacked2/cluster_qsub.sh \
--use-singularity \
--singularity-args "-B /SAN/vyplab/:/SAN/vyplab/" \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} {cluster.submission_string}" \
-j 40 \
--nolock \
--rerun-incomplete \
--latency-wait 100