import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "../rules/helpers.py"

samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name']))

print(SAMPLE_NAMES)

SPECIES = config["species"]
GTF = config['gtf']

BASES, CONTRASTS = return_bases_and_contrasts()




include: "../rules/create_cryptic_exons_majiq_deltas.smk"


rule allParseFinalBeds:
    input:
        expand(os.path.join(OUTPUT_CRYPTIC_EXONS,"{bse}_{contrast}_cryptic_exons.bed"),zip, bse = BASES,contrast = CONTRASTS)
