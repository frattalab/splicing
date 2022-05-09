import pandas as pd
import os
import subprocess
import yaml

configfile: "config/config.yaml"

include: "helpers.py"

#reading in the samples and dropping the samples to be excluded in order to get a list of sample names

samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
SAMPLE_NAMES_NOPERIODS = list(set(samples2['sample_name']))
SJ_NAMES  = list(set(samples2['sample_name']))
print(SAMPLE_NAMES)
GROUPS = list(set(samples2['group']))

BASES, CONTRASTS = return_bases_and_contrasts()
MAJIQ_DIR = get_output_dir(config['project_top_level'], config['majiq_top_level'])

rule allParse:
    input:
        expand(os.path.join(MAJIQ_DIR,"modulizers","{bse}-{contrast}", "junctions.tsv"),zip, bse = BASES,contrast = CONTRASTS),
        expand(os.path.join(MAJIQ_DIR,"psi_voila_tsv_single",'{sample}' + "_parsed.csv"), sample = SAMPLE_NAMES),
        expand(os.path.join(MAJIQ_DIR,"psi_voila_tsv",'{group}' + "_parsed.csv"), group = GROUPS),
        expand(os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv","{bse}-{contrast}" + "_parsed_psi.tsv"),zip, bse = BASES,contrast = CONTRASTS),
        

include: "../rules/majiq_modulizer.smk"
include: "../rules/majiq_parse_deltas.smk"
include: "../rules/majiq_annotate_junctions.smk"
