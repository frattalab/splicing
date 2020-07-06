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


rule all:
    input:
        config['majiq_top_level'] + config['run_name'] + "_majiqConfig.tsv",
        expand(os.path.join(config['majiq_top_level'],"builder",'{name}' + ".majiq"),name = SAMPLE_NAMES),

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
        os.path.join(config['majiq_top_level'],"builder/splicegraph.sql"),
        os.path.join(config['majiq_top_level'],"builder/builder_done")
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
        touch {params.majiq_builder_output}/builder_done
        """
