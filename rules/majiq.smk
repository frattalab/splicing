import pandas as pd
import os
import subprocess
import yaml
# note to self - probaly should not be hard coding threads like i do
configfile: "config/config.yaml"
include: "helpers.py"
localrules: create_majiq_config_file
#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))

BASES, CONTRASTS = return_bases_and_contrasts()
print(BASES)
print(CONTRASTS)

rule all:
    input:
        expand(config['majiq_top_level'] + "delta_psi/" + "{bse}_{contrast}" + ".deltapsi.tsv",zip, bse = BASES,contrast = CONTRASTS),
        expand(os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}" + ".psi.tsv"),zip, bse = BASES,contrast = CONTRASTS),
        config['majiq_top_level'] + config['run_name'] + "_majiqConfig.tsv",
        expand(os.path.join(config['majiq_top_level'],"builder",'{name}' + ".majiq"),name = SAMPLE_NAMES),
        expand(os.path.join(config['majiq_top_level'],"psi",'{group}' + ".psi.voila"),group = GROUPS)

# # this rule creats the majiq configuration file that is required uses a helper function
rule create_majiq_config_file:
    input:
        config['sample_csv_path']
    output:
        majiq_config = config['majiq_top_level'] + config['run_name'] + "_majiqConfig.tsv"
    # use config.yaml and samples.tsv to make MAJIQ file
    run:
        conditions_bams_parsed = parse_sample_csv_majiq(config['sample_csv_path'])
        options = [
            "[info]",
            "readlen=" + str(config['read_len']),
            "bamdirs=" + config['bam_dir'],
            "genome=" + config['genome_refcode'],
            "strandness=" + config['strand_code'],
            "[experiments]"]

        options += conditions_bams_parsed
        with open(output.majiq_config,"w") as outfile:
            for opt in options:
                outfile.write(opt + "\n")

rule majiq_build:
    input:
        majiq_config_file = config['majiq_top_level'] + config['run_name'] + "_majiqConfig.tsv"
    output:
        expand(os.path.join(config['majiq_top_level'],"builder",'{name}' + ".majiq"),name = SAMPLE_NAMES),
        os.path.join(config['majiq_top_level'],"builder/splicegraph.sql")
    threads:
            4
    params:
        majiq_path = config['majiq_path'],
        gff3 = config['gff3'],
        majiq_builder_output = os.path.join(config['majiq_top_level'],"builder"),
        majiq_extra_parameters = return_parsed_extra_params(config['majiq_extra_parameters'])
    shell:
        """
        mkdir -p {params.majiq_builder_output}
        {params.majiq_path} build {params.gff3} -c {input.majiq_config_file} -j {threads} -o {params.majiq_builder_output}{params.majiq_extra_parameters}
        """

rule majiq_psi:
    input:
    #this is always calling from the column named 'group' in the sample csv file
        group_majiq = lambda wildcards: majiq_files_by_group(wildcards.group)
    output:
        voila = os.path.join(config['majiq_top_level'],"psi",'{group}' + ".psi.voila"),
        tsv = os.path.join(config['majiq_top_level'],"psi",'{group}' + ".psi.tsv")
    params:
        majiq_path = config['majiq_path'],
        psi_output_folder = os.path.join(config['majiq_top_level'],"psi"),
    threads:
        16
    shell:
        """
        mkdir -p {params.psi_output_folder}
        {params.majiq_path} psi {input.group_majiq} -j {threads} -o {params.psi_output_folder} -n {wildcards.group}
        """
rule majiq_delta_psi:
    input:
        majiq_config_file = config['majiq_top_level'] + config['run_name'] + "_majiqConfig.tsv",
        finished_splicegraph = os.path.join(config['majiq_top_level'],"builder/splicegraph.sql"),
        base_group_majiq = lambda wildcards: majiq_files_from_contrast(wildcards.bse),
        contrast_group_majiq = lambda wildcards: majiq_files_from_contrast(wildcards.contrast)
    output:
        os.path.join(config['majiq_top_level'],"delta_psi","{bse}_{contrast}" + ".deltapsi.tsv"),
        os.path.join(config['majiq_top_level'],"delta_psi","{bse}_{contrast}" + ".deltapsi.voila")
    params:
        majiq_path = config['majiq_path'],
        delta_psi_output_folder = os.path.join(config['majiq_top_level'],"delta_psi"),
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
        voila_file = lambda wildcards: os.path.join(config['majiq_top_level'],"delta_psi","{bse}_{contrast}" + ".deltapsi.voila")
    output:
        tsv = os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}" + ".psi.tsv")
    params:
        voila_path = config['voila_path'],
        psi_output_folder = os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv"),
        splice_graph = os.path.join(config['majiq_top_level'],"builder", "splicegraph.sql"),
        extra_voila_paramters = return_parsed_extra_params(config['extra_voila_parameters'])
    shell:
        """
        mkdir -p {params.psi_output_folder}
        {params.voila_path} tsv {params.splice_graph} {input.voila_file} -f {output.tsv} {params.extra_voila_paramters}
        """
