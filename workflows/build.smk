import pandas as pd
import os
import subprocess
import yaml
# note to self - probaly should not be hard coding threads like i do
configfile: "config/config.yaml"
include: "../rules/helpers.py"
localrules: create_majiq_config_file
#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))


rule all:
    input:
        config['majiq_top_level'] + config['run_name'] + "_majiqConfig.tsv",
        expand(os.path.join(config['majiq_top_level'],"builder",'{name}' + ".majiq"),name = SAMPLE_NAMES)

include: "../rules/majiq_build.smk"
