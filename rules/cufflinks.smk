import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"

rule merged_cuff:
    input:
        cuff_link_folder + "/merged_asm.gtf"

rule cufflinks:
    input:
        "{sample}.bam"
    output:
        cuff_link_folder + sample.gtf
    shell:
    """
    cufflinks ./test.ctrl2.bam --library-type fr-firststrand --GTF-guide /SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v31.no_chr.annotation.gtf -o {}
    """

rule cuffmerge:
    input:
        expand(cuff_link_folder + "{samples}.gtf")
    output:
        cuff_link_folder
    shell:
    """
    """
