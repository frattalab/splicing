

import os
project_dir = "/SAN/vyplab/alb_projects/alb_projects/data/ibm_geo/"

bam_dir = os.path.join(project_dir,"STAR_aligned/")

output_dir = os.path.join(project_dir,"picardmetrics/")

SAMPLES, = glob_wildcards(bam_dir + "{sample}.Aligned.sorted.out.bam")

rule collect_insert_size:
    input:
        sample_bam = lambda wildcards: bam_dir + "{sample}.Aligned.sorted.out.bam"
    output:
    # SBMA_2.Aligned.sorted.out.bam.insert_size_metrics.txt
        insertmetrics = os.path.join(project_dir,"picardmetrics/")
    shell:
        """
        mkdir -p {output_dir}
        java -jar picard.jar CollectInsertSizeMetrics \
        I=input.bam \
        O= {output.insertmetrics} \
        H=insert_size_histogram.pdf
        """
