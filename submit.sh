#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=8G,h_rt=20:00:00,tmem=8G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y

WORKFLOW="workflows/${1}.smk"

if [ "$2" != "" ]; then
    RUN_NAME=$1_$2
else
    RUN_NAME=$1
fi

FOLDER=submissions/$(date +"%Y%m%d%H%M")

mkdir -p ${FOLDER}
cp config/config.yaml ${FOLDER}/${RUN_NAME}_config.yaml
cp config/comparisons.yaml ${FOLDER}/${RUN_NAME}_comparisons.yaml


snakemake -s ${WORKFLOW} \
--use-conda \
--conda-prefix "/SAN/vyplab/vyplab_reference_genomes/conda_envs/" \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o $FOLDER {cluster.submission_string}" \
-j 40 \
--nolock \
--rerun-incomplete \
--latency-wait 100 \
