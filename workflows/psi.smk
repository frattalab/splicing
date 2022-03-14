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



rule allPSIwf:
    input:
        expand(os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv","{bse}-{contrast}" + ".psi.tsv"),zip, bse = BASES,contrast = CONTRASTS),
        expand(os.path.join(MAJIQ_DIR,"psi",'{group}' + ".psi.voila"),group = GROUPS),
        expand(os.path.join(MAJIQ_DIR,"psi_voila_tsv",'{group}' + ".psi.tsv"),group = GROUPS),
        expand(os.path.join(MAJIQ_DIR,"psi_single","{sample}" + config['bam_suffix'] + ".psi.voila"), sample = SAMPLE_NAMES_NOPERIODS),
        expand(os.path.join(MAJIQ_DIR,"psi_voila_tsv_single",'{sample}' + config['bam_suffix'] + ".psi.tsv"), sample = SAMPLE_NAMES_NOPERIODS)


include: "../rules/majiq_psi.smk"
include: "../rules/majiq_parse_deltas.smk"
include: "../rules/majiq_annotate_junctions.smk"
