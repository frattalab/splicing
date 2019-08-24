#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1G,h_rt=48:00:00,tmem=1G
#$ -pe smp 16

# join stdout and stderr output
#$ -j y
#$ -R y
#$ -N majic_no_snake
export LD_LIBRARY_PATH=/share/apps/genomics/htslib-1.9/lib:$LD_LIBRARY_PATH
source /share/apps/examples/source_files/python/python-3.6.9.source

/share/apps/python-3.6.9/bin/majiq build /SAN/vyplab/vyplab_reference_genomes/annotation/mouse/ensembl.mm10.no_chr.gff3 \
-c /SAN/vyplab/alb_projects/data/4su_tdp_f210i/majiq/4su_mef_f210_wt.tsv -j 16 \
-o /SAN/vyplab/alb_projects/data/4su_tdp_f210i/majiq/builder_simple --simplify
