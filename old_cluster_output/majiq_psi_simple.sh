#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1G,h_rt=48:00:00,tmem=1G
#$ -pe smp 16

# join stdout and stderr output
#$ -j y
#$ -R y
#$ -N majiq_psi_simple
export LD_LIBRARY_PATH=/share/apps/genomics/htslib-1.9/lib:$LD_LIBRARY_PATH
source /share/apps/examples/source_files/python/python-3.6.9.source

/share/apps/python-3.6.9/bin/majiq psi /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/builder_simple/Ctrl*.majiq -j 16 -o /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/psi_output_simple -n control

# /share/apps/python-3.6.9/bin/majiq psi /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/builder_simple/TDP_*.majiq -j 16 -o /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/psi_output_simple -n tdp
#
# /share/apps/python-3.6.9/bin/majiq psi /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/builder_simple/C9_*.majiq -j 16 -o /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/psi_output_simple -n c9
#
# /share/apps/python-3.6.9/bin/majiq psi /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/builder_simple/Fus_*.majiq -j 16 -o /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/psi_output_simple -n fus
#
# /share/apps/python-3.6.9/bin/majiq psi /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/builder_simple/SBMA_*.majiq -j 16 -o /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/psi_output_simple -n sbma
#
# /share/apps/python-3.6.9/bin/majiq psi /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/builder_simple/SOD_*.majiq -j 16 -o /SAN/vyplab/alb_projects/data/muscle/analysis/majiq/psi_output_simple -n sod
