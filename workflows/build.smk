import pandas as pd
import os
import subprocess
import yaml
# note to self - probaly should not be hard coding threads like i do
configfile: "config/config.yaml"
include: "../rules/helpers.py"
localrules: create_majiq_config_file
#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))
MAJIQ_DIR = get_output_dir(config['project_top_level'], config['majiq_top_level'])


rule allBuild:
    input:
        MAJIQ_DIR + config['run_name'] + "_majiqConfig.tsv",
        expand(os.path.join(MAJIQ_DIR,"builder",'{name}' + ".sj"),name = SAMPLE_NAMES),
        os.path.join(MAJIQ_DIR,"splicegraph.sql")


include: "../rules/majiq_build.smk"
