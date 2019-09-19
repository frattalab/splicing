#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1G,h_rt=48:00:00,tmem=1G
#$ -pe smp 8

# join stdout and stderr output
#$ -j y
#$ -R y
#$ -N cuff_link_try_again

# cufflinks /home/annbrown/muscle/cuffdiff_gtf/Ctrl1.mvk.sam --library-type fr-firststrand --GTF-guide /SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v31.no_chr.annotation.gtf -o /home/annbrown/muscle/cuffdiff_gtf/Ctrl1
#
# cufflinks /home/annbrown/muscle/cuffdiff_gtf/Fus1.mvk.sam --library-type fr-firststrand --GTF-guide /SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v31.no_chr.annotation.gtf -o /home/annbrown/muscle/cuffdiff_gtf/Fus1

cuffmerge /home/annbrown/pipelines/splicing/gtf_list.txt --ref-gtf /SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v31.no_chr.annotation.gtf --keep-tmp
-p 8
