#$ -S /bin/bash
#$ -l h_vmem=4G,tmem=4G
#$ -l h_rt=72:00:00
#$ -pe smp 2
#$ -R y
#$ -j y

source activate sgseq
Rscript /SAN/vyplab/alb_projects/pipelines/splicing/scripts/sg_seq_create_anno.R
