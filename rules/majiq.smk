import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"
include: "helpers.py"

#reading in the samples and dropping the sampels to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))


rule top:
    input:
        config['majiq_top_level'] + config['run_name'] + "_majiqConfig.tsv",
        os.path.join(config['majiq_top_level'],"builder","majiq.log"),
        expand(os.path.join(config['majiq_top_level'],"psi",'{group}' + ".psi.voila"),group = GROUPS)

rule create_majiq_config_file:
    input:
        config['sample_csv_path']
    output:
        majiq_config = config['majiq_top_level'] + config['run_name'] + "_majiqConfig.tsv"
    params:
        conditions_bams_parsed = parse_sample_csv_majiq(config['sample_csv_path'])
    # use config.yaml and samples.tsv to make MAJIQ file
    run:
        options = [
            "[info]",
            "readlen=" + str(config['read_len']),
            "bamdir=" + config['bam_dir'],
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
