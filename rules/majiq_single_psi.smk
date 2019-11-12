import pandas as pd
import os
import subprocess
import yaml

configfile: "config/config.yaml"
include: "helpers.py"

#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))

rule top:
    input:
        expand(os.path.join(config['majiq_top_level'],"psi_single",'{sample}' + ".psi.voila"),sample = SAMPLE_NAMES),
        expand(os.path.join(config['majiq_top_level'],"psi_voila_tsv_single",'{sample}' + ".psi.tsv"), sample = SAMPLE_NAMES)

rule majiq_psi:
    input:
        group_majiq = lambda wildcards: os.path.join(config['majiq_top_level'],"builder",wildcards.sample + ".majiq")
    output:
        voila = os.path.join(config['majiq_top_level'],"psi_single",'{sample}' + ".psi.voila"),
        tsv = os.path.join(config['majiq_top_level'],"psi_single",'{sample}' + ".psi.tsv")
    params:
        majiq_path = config['majiq_path'],
        psi_output_folder = os.path.join(config['majiq_top_level'],"psi_single"),
    threads:
        4
    shell:
        """
        mkdir -p {params.psi_output_folder}
        {params.majiq_path} psi {input.group_majiq} -j {threads} -o {params.psi_output_folder} -n {wildcards.sample}
        """
rule majiq_psi_tsv:
    input:
    #this is always calling from the column named 'group' in the sample csv file
        voila_file = lambda wildcards: os.path.join(config['majiq_top_level'],"psi_single",'{sample}' + ".psi.voila")
    output:
        tsv = os.path.join(config['majiq_top_level'],"psi_voila_tsv_single",'{sample}' + ".psi.tsv")
    params:
        voila_path = config['voila_path'],
        psi_output_folder = os.path.join(config['majiq_top_level'],"psi_voila_tsv_single"),
        splice_graph = os.path.join(config['majiq_top_level'],"builder", "splicegraph.sql")
    shell:
        """
        mkdir -p {params.psi_output_folder}
        {params.voila_path} tsv {params.splice_graph} {input.voila_file} -f {output.tsv}
        """
