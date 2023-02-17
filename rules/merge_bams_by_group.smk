import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "../rules/helpers.py"
##############################
##### STOLEN FROM https://github.com/bioinformatics-core-shared-training/RNAseq_March_2019/tree/master/
##### Final output is a merged gtf of all samples containing transcripts that
##### were found in the samples but NOT the input GTF
##############################
#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name']))
BASES, CONTRASTS = return_bases_and_contrasts()
ALL_GROUPS = list(set(BASES,CONTRASTS))

print(ALL_GROUPS)

SPECIES = config["species"]
GTF = config['gtf']

#make sure the output folder for STAR exists before running anything
bam_dir = get_output_dir(config["project_top_level"], config['bam_dir'])
bam_dir_temporary = get_output_dir(config["project_top_level"], 'merged_bams')

stringtie_outdir = get_output_dir(config["project_top_level"], config['stringtie_output'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])

# RULE ORDER DIRECTIVE
# if paired end, use the paired end rule to run, if single end use the single end rule to run
ruleorder: merge_bam_base > merge_bam_contrast


rule allMerging:
    input:
        expand(os.path.join(bam_dir_temporary,"{bse}.bam"), bse = BASES),
        expand(os.path.join(bam_dir_temporary,"{contrast}.bam"), contrast = CONTRASTS)


rule merge_bam_base:
    input:
        base_group_stringtie = lambda wildcards: file_path_list(wildcards.bse,bam_dir,config['bam_suffix'] + '.bam')
    output:
        bam= os.path.join(bam_dir_temporary,"{bse}.bam"),
        bai= os.path.join(bam_dir_temporary,"{bse}.bam.bai")
    threads: 10
    conda:
        "../envs/stringtie.yml"
    shadow:
        "shallow"
    shell:
        "samtools merge {output.bam} {input.base_group_stringtie} --threads {threads};samtools index {output.bam} {output.bai} 2> {log} "

rule merge_bam_contrast:
    input:
        contrast_group_stringtie = lambda wildcards: file_path_list(wildcards.contrast,bam_dir,config['bam_suffix'] + '.bam')
    output:
        bam= os.path.join(bam_dir_temporary,"{contrast}.bam"),
        bai= os.path.join(bam_dir_temporary,"{contrast}.bam.bai")
    threads: 10
    conda:
        "../envs/stringtie.yml"
    shadow:
        "shallow"
    shell:
        "samtools merge {output.bam} {input.contrast_group_stringtie} --threads {threads};samtools index {output.bam} {output.bai} 2> {log} "
