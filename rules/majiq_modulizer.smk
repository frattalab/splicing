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
SAMPLE_NAMES_NOPERIODS = list(set(samples2['sample_name']))

GROUPS = list(set(samples2['group']))

BASES, CONTRASTS = return_bases_and_contrasts()
print(BASES)
print(CONTRASTS)
MAJIQ_DIR = get_output_dir(config['project_top_level'], config['majiq_top_level'])

rule allPSI:
    input:
        expand(os.path.join(MAJIQ_DIR,"modulizers","{bse}-{contrast}", "junctions.tsv"),zip, bse = BASES,contrast = CONTRASTS)

rule majiq_delta_modulizer:
    input:
        os.path.join(MAJIQ_DIR,"delta_psi","{bse}-{contrast}" + ".deltapsi.voila")
    output:
        os.path.join(MAJIQ_DIR,"modulizers","{bse}-{contrast}", "junctions.tsv"),
    params:
        voila_path = config['voila_path'],
        delta_mod_output_folder = os.path.join(MAJIQ_DIR,"modulizers","{bse}-{contrast}"),
        splice_graph = os.path.join(MAJIQ_DIR,"builder", "splicegraph.sql")
    threads:
        8
    shell:
        """
        mkdir -p {params.delta_psi_output_folder}
        {params.voila_path} modulize {params.splice_graph} \
        {input} \
        -d {params.delta_mod_output_folder} -j {threads}
        """
