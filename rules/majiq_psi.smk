import pandas as pd
import os
import subprocess
import yaml
# note to self - probaly should not be hard coding threads like i do
configfile: "config/config.yaml"
include: "helpers.py"
localrules: create_majiq_config_file
#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))

BASES, CONTRASTS = return_bases_and_contrasts()
print(BASES)
print(CONTRASTS)
MAJIQ_DIR = get_output_dir(config['project_top_level'], config['majiq_top_level'])

rule allPSI:
    input:
        expand(os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv","{bse}-{contrast}" + ".psi.tsv"),zip, bse = BASES,contrast = CONTRASTS),
        expand(os.path.join(MAJIQ_DIR,"psi",'{group}' + ".psi.voila"),group = GROUPS),
        # expand(os.path.join(MAJIQ_DIR,"psi_single",'{sample}' + ".psi.voila"),sample = SAMPLE_NAMES),
        # expand(os.path.join(MAJIQ_DIR,"psi_voila_tsv_single",'{sample}' + ".psi.tsv"), sample = SAMPLE_NAMES)

rule majiq_psi:
    input:
    #this is always calling from the column named 'group' in the sample csv file
        group_majiq = lambda wildcards: majiq_files_by_group(wildcards.group)
    output:
        voila = os.path.join(MAJIQ_DIR,"psi",'{group}' + ".psi.voila"),
        tsv = os.path.join(MAJIQ_DIR,"psi",'{group}' + ".psi.tsv")
    params:
        majiq_path = config['majiq_path'],
        psi_output_folder = os.path.join(MAJIQ_DIR,"psi"),
    threads:
        16
    shell:
        """
        mkdir -p {params.psi_output_folder}
        {params.majiq_path} psi {input.group_majiq} -j {threads} -o {params.psi_output_folder} -n {wildcards.group}
        """
rule majiq_delta_psi:
    input:
        majiq_config_file = MAJIQ_DIR + config['run_name'] + "_majiqConfig.tsv",
        base_group_majiq = lambda wildcards: majiq_files_from_contrast(wildcards.bse),
        contrast_group_majiq = lambda wildcards: majiq_files_from_contrast(wildcards.contrast)
    output:
        os.path.join(MAJIQ_DIR,"delta_psi","{bse}-{contrast}" + ".deltapsi.tsv"),
        os.path.join(MAJIQ_DIR,"delta_psi","{bse}-{contrast}" + ".deltapsi.voila")
    params:
        majiq_path = config['majiq_path'],
        delta_psi_output_folder = os.path.join(MAJIQ_DIR,"delta_psi"),
        majiq_psi_extra_parameters = return_parsed_extra_params(config['majiq_psi_extra_parameters'])
    threads:
        8
    shell:
        """
        mkdir -p {params.delta_psi_output_folder}
        {params.majiq_path} deltapsi -grp1 {input.base_group_majiq} -grp2 {input.contrast_group_majiq} -j {threads} -o {params.delta_psi_output_folder}{params.majiq_psi_extra_parameters} --name {wildcards.bse} {wildcards.contrast}
        """
rule majiq_delta_psi_tsv:
    input:
    #this is always calling from the column named 'group' in the sample csv file
        voila_file = lambda wildcards: os.path.join(MAJIQ_DIR,"delta_psi","{bse}-{contrast}" + ".deltapsi.voila")
    output:
        tsv = os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv","{bse}-{contrast}" + ".psi.tsv")
    params:
        voila_path = config['voila_path_old'],
        psi_output_folder = os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv"),
        splice_graph = os.path.join(MAJIQ_DIR,"builder", "splicegraph.sql"),
        extra_voila_paramters = return_parsed_extra_params(config['extra_voila_parameters'])
    shell:
        """
        mkdir -p {params.psi_output_folder}
        {params.voila_path} tsv {params.splice_graph} {input.voila_file} -f {output.tsv} {params.extra_voila_paramters}
        """

# rule majiq_single_psi:
#     input:
#         group_majiq = lambda wildcards: os.path.join(MAJIQ_DIR,"builder",wildcards.sample + ".majiq")
#     output:
#         voila = os.path.join(MAJIQ_DIR,"psi_single",'{sample}.replace(config['bam_suffix'],"") + ".psi.voila"),
#         tsv = os.path.join(MAJIQ_DIR,"psi_single",'{sample}' + ".psi.tsv")
#     params:
#         majiq_path = config['majiq_path'],
#         psi_output_folder = os.path.join(MAJIQ_DIR,"psi_single"),
#         test = lambda wildcards: wildcards.sample.replace(config['bam_suffix'],"")
#     threads:
#         4
#     shell:
#         """
#         mkdir -p {params.psi_output_folder}
#         echo {params.test}
#         {params.majiq_path} psi {input.group_majiq} -j {threads} -o {params.psi_output_folder} -n {params.test}
#       """

rule majiq_psi_tsv:
    input:
    #this is always calling from the column named 'group' in the sample csv file
        voila_file = lambda wildcards: os.path.join(MAJIQ_DIR,"psi_single",'{sample}' + ".psi.voila")
    output:
        tsv = os.path.join(MAJIQ_DIR,"psi_voila_tsv_single",'{sample}' + ".psi.tsv")
    params:
        voila_path = config['voila_path_old'],
        psi_output_folder = os.path.join(MAJIQ_DIR,"psi_voila_tsv_single"),
        splice_graph = os.path.join(MAJIQ_DIR,"builder", "splicegraph.sql")
    shell:
        """
        mkdir -p {params.psi_output_folder}
        {params.voila_path} tsv {params.splice_graph} {input.voila_file} -f {output.tsv}
        """
