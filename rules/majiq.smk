import pandas as pd
import os
import subprocess
import yaml

configfile: "config/config.yaml"
include: "helpers.py"

#reading in the samples and dropping the sampels to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))
BASES, CONTRASTS = return_bases_and_contrasts()

rule top:
    input:
        config['majiq_top_level'] + config['run_name'] + "_majiqConfig.tsv",
        os.path.join(config['majiq_top_level'],"builder","majiq.log"),
        expand(os.path.join(config['majiq_top_level'],"psi",'{group}' + ".psi.voila"),group = GROUPS),
        expand(os.path.join(config['majiq_top_level'],"psi",'{base}_{contrast}' + '.tsv'),zip,base = BASES,contrast = CONTRASTS)

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
        os.path.join(config['majiq_top_level'],"builder","majiq.log"),
        expand(os.path.join(config['majiq_top_level'],"builder",'{name}' + ".majiq"),name = SAMPLE_NAMES)
    threads:
            16
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
        base_group_majiq = lambda wildcards: majiq_files_from_contrast(wildcards.base),
        contrast_group_majiq = lambda wildcards: majiq_files_from_contrast(wildcards.contrast)
    output:
        os.path.join(config['majiq_top_level'],"delta_psi",'{base}_{contrast}' + '.tsv')
    params:
        majiq_path = config['majiq_path'],
        delta_psi_output_folder = os.path.join(config['majiq_top_level'],"delta_psi"),
        wildcards.
    threads:
        16
    shell:
        """
        mkdir -p {params.delta_psi_output_folder}
        {params.majiq_path} deltapsi -grp1 {input.control_group_majiq} grp2 {input.control_group_majiq} -j {threads} -o {params.delta_psi_output_folder} --name {wildcards.base} {wildcards.contrast}
        """
