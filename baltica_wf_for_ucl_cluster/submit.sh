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
CONFIGFILE="/SAN/vyplab/first_weeks/BDNF_4SU_FULL/analysis/baltica/config.yml"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

##Suggested run order
#junctionseq.smk  
#leafcutter.smk  
##You HAVE to do build before PSI
#majiq_build.smk  
#majiq_psi.smk
#rmats.smk
#stringtie.smk
##analysis is the last thing which is run

#analysis.smk  

snakemake -s junctionseq.smk \
--configfile $CONFIGFILE \
--cluster-config $SCRIPT_DIR/cluster.yaml \
--jobscript $SCRIPT_DIR/cluster_qsub.sh \
--use-singularity \
--singularity-args "-B /SAN/vyplab/:/SAN/vyplab/" \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} {cluster.submission_string}" \
-j 40 \
--nolock \
--rerun-incomplete \
--latency-wait 100
