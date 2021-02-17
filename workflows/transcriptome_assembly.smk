import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "../rules/helpers.py"
localrules: compose_gtf_list
##############################
##### STOLEN FROM https://github.com/bioinformatics-core-shared-training/RNAseq_March_2019/tree/master/
##### Final output is a merged gtf of all samples containing transcripts that
##### were found in the samples but NOT the input GTF
##############################
#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name']))

print(SAMPLE_NAMES)

SPECIES = config["species"]
GTF = config['gtf']

#make sure the output folder for STAR exists before running anything
bam_dir = get_output_dir(config["project_top_level"], config['bam_dir'])
stringtie_outdir = get_output_dir(config["project_top_level"], config['stringtie_output'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])

rule assemble_both:
    input:
        expand(scallop_outdir + '{sample}' + ".gtf", sample = SAMPLE_NAMES),
        expand(scallop_outdir + "gffall.{sample}.gtf.tmap",sample = SAMPLE_NAMES),
        expand(scallop_outdir + "{sample}.unique.gtf",sample = SAMPLE_NAMES),
        os.path.join(scallop_outdir,"scallop_merged.gtf"),
        expand(stringtie_outdir + "{sample}.assemble.gtf", sample = SAMPLE_NAMES),
        expand(stringtie_outdir + "gffall.{sample}.gtf.tmap",sample = SAMPLE_NAMES),
        expand(stringtie_outdir + "{sample}.unique.gtf",sample = SAMPLE_NAMES),
        os.path.join(stringtie_outdir,"stringtie_merged.gtf")


include: "../rules/scallop_sample.smk"
include: "../rules/stringtie_sample.smk"
