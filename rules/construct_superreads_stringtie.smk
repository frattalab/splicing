import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "helpers.py"
##############################
##### Final output is a merged gtf of all samples containing transcripts that
##### were found in the samples but NOT the input GTF
##############################


GTF = config['gtf']

#make sure the output folder for STAR exists before running anything
bam_dir = get_output_dir(config["project_top_level"], config['bam_dir'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])



rule index:
    input:
        WORKING_DIR + "genome/genome.fasta"
    output:
        [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    conda:
        "../envs/stringtie_superreads.yaml"
    message:
        "indexing genome"
    params:
        WORKING_DIR + "genome/genome"
    threads: 10
    shell:
        "hisat2-build -p {threads} {input} {params} --quiet"
build the gmap indexes
