#! /bin/bash
# this file is snakemake.sh
module load snakemake  || exit 1
module load majiq
module load R

WORKFLOW="workflows/${1}.smk"

if [ "$2" != "" ]; then
    RUN_NAME=$1_$2
else
    RUN_NAME=$1
fi

FOLDER=submissions/$(date +"%Y%m%d%H%M")

mkdir -p ${FOLDER}
cp config/config.yaml ${FOLDER}/${RUN_NAME}_config.yaml

snakemake -s ${WORKFLOW} \
-pr \
--jobs 40 \
--keep-going \
--cluster-config config/cluster_slurm.yaml \
--cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus}" \
--nolock \
--rerun-incomplete \
--latency-wait 100 \
--use-conda
