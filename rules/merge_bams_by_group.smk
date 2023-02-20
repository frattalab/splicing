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
ALLGROUP = list(set(BASES + CONTRASTS))
print(ALLGROUP)


SPECIES = config["species"]
GTF = config['gtf']

#make sure the output folder for STAR exists before running anything
bam_dir = get_output_dir(config["project_top_level"], config['bam_dir'])
bam_dir_temporary = get_output_dir(config["project_top_level"], 'merged_bams')

stringtie_outdir = get_output_dir(config["project_top_level"], config['stringtie_output'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])



rule allMerging:
    input:
        expand(os.path.join(bam_dir_temporary,"{grp}.bam"), grp = ALLGROUP)


rule merge_bam_groups:
    input:
        base_group_stringtie = lambda wildcards: file_path_list(wildcards.grp,bam_dir,config['bam_suffix'] + '.bam')
    output:
        bam= os.path.join(bam_dir_temporary,"{grp}.bam"),
        bai= os.path.join(bam_dir_temporary,"{grp}.bam.bai")
    threads: 10
    shadow:
        "shallow"
    shell:
        "samtools merge {output.bam} {input.base_group_stringtie} --threads {threads};samtools index {output.bam} {output.bai} 2> {log} "

