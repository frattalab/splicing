import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "helpers.py"
localrules: compose_gtf_list

#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))
print(SAMPLE_NAMES)

rule all_scallop:
    input:
        expand(os.path.join(config['top_level_project_folder'],"scallop_output/",'{sample}' + ".gtf"), sample = SAMPLE_NAMES),
        os.path.join(config['top_level_project_folder'],"scallop_output/","scallop_merged.gtf"),
        os.path.join(config['top_level_project_folder'],"scallop_output/","gffall.scallop_merged.gtf.map")


rule scallop_per_samp:
    input:
        bam_file = lambda wildcards: config['bam_dir'] + '{sample}' + ".bam"
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    output:
        os.path.join(config['top_level_project_folder'],"scallop_output/",'{sample}' + ".gtf")
    params:
        scallop_path = config['scallop_path'],
        scallop_out_folder = os.path.join(config['top_level_project_folder'],"scallop_output/"),
        scallop_extra_config = return_parsed_extra_params(config['scallop_extra_parameters'])
    shell:
        """
        mkdir -p {params.scallop_out_folder}
        {params.scallop_path} -i {input.bam_file} -o {output} {params.scallop_extra_config}
        """

rule compose_gtf_list:
    input:
        expand(os.path.join(config['top_level_project_folder'],"scallop_output/",'{sample}.gtf'), sample=SAMPLE_NAMES)
    output:
        txt = os.path.join(config['top_level_project_folder'],"scallop_output/","gtf_list.txt")
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)
rule merge_scallop_gtfs:
    input:
        gtf_list = os.path.join(config['top_level_project_folder'],"scallop_output/","gtf_list.txt")
    output:
        merged_gtf = os.path.join(config['top_level_project_folder'],"scallop_output/","scallop_merged.gtf")
    params:
        gtfmerge = '/SAN/vyplab/alb_projects/tools/rnaseqtools-1.0.3/gtfmerge/gtfmerge'
    shell:
        """
        {params.gtfmerge} union {input.gtf_list} {output.merged_gtf} -t 2 -n
        """
rule compare_reference:
    input:
        merged_gtf = os.path.join(config['top_level_project_folder'],"scallop_output/","scallop_merged.gtf")
    output:
        os.path.join(config['top_level_project_folder'],"scallop_output/","gffall.scallop_merged.gtf.tmap")
    params:
        ref_gtf = config['gtf'],
        gffcompare = "/SAN/vyplab/alb_projects/tools/gffcompare-0.11.6.Linux_x86_64/gffcompare"
    shell:
        """
        {params.gffcompare} -o gffall -r {params.ref_gtf} {input.merged_gtf}
        """
