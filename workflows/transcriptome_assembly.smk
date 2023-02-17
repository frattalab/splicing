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
BASES, CONTRASTS = return_bases_and_contrasts()
print(SAMPLE_NAMES)

SPECIES = config["species"]
GTF = config['gtf']

#make sure the output folder for STAR exists before running anything
bam_dir = get_output_dir(config["project_top_level"], config['bam_dir'])
stringtie_outdir = get_output_dir(config["project_top_level"], config['stringtie_output'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])


include: "../rules/merge_bams_by_group.smk"
include: "../rules/scallop_sample.smk"
include: "../rules/stringtie_sample.smk"

rule assemble_all:
    input:
        expand(scallop_outdir + '{sample}' + ".gtf", sample = SAMPLE_NAMES),
        os.path.join(scallop_outdir, "scallop_merged.unique.gtf"),
        os.path.join(scallop_outdir, "gffall.scallop_merged.gtf.tmap"),
        expand(os.path.join(scallop_outdir,"{bse}.scallop_merged.gtf"), bse = BASES),
        expand(os.path.join(scallop_outdir,"{contrast}.scallop_merged.gtf"), contrast = CONTRASTS),
        expand(stringtie_outdir + "{sample}.assemble.gtf", sample = SAMPLE_NAMES),
        os.path.join(stringtie_outdir, "stringtie_merged.unique.gtf"),
        os.path.join(stringtie_outdir,"stringtie_merged.gtf"),
        expand(os.path.join(scallop_outdir,"{bse}.scallop_merged.gtf"), bse = BASES),
        expand(os.path.join(scallop_outdir,"{contrast}.scallop_merged.gtf"), contrast = CONTRASTS),
        os.path.join(config["project_top_level"],"all_assemblers_merged.gtf"),
        expand(os.path.join(stringtie_outdir,"{bse}.stringtie_merged.gtf"), bse = BASES),
        expand(os.path.join(stringtie_outdir,"{contrast}.stringtie_merged.gtf"), contrast = CONTRASTS)

rule compose_gtf_list_all_assemblers:
    input:
        os.path.join(stringtie_outdir,"stringtie_merged.gtf"),
        os.path.join(scallop_outdir, "scallop_merged.gtf"),
    output:
        txt = os.path.join(config["project_top_level"],"gtf_list_all_assemblers.txt")
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)

rule merge_all_assemblies:
    input:
        gtf_list = os.path.join(config["project_top_level"],"gtf_list_all_assemblers.txt")
    output:
        merged_gtf = os.path.join(config["project_top_level"],"all_assemblers_merged.gtf")
    params:
        gtfmerge = config['gtfmerge']
    shell:
        """
        {params.gtfmerge} union {input.gtf_list} {output.merged_gtf} -t 2 -n
        """
