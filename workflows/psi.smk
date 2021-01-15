import pandas as pd
import os
import subprocess
import yaml
# note to self - probaly should not be hard coding threads like i do
configfile: "config/config.yaml"
include: "../rules/helpers.py"
#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))

BASES, CONTRASTS = return_bases_and_contrasts()
print(BASES)
print(CONTRASTS)
import pandas as pd
import os
import subprocess
import yaml


rule allPSIwf:
    input:
        expand(os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}" + ".psi.tsv"),zip, bse = BASES,contrast = CONTRASTS),
        expand(os.path.join(config['majiq_top_level'],"psi",'{group}' + ".psi.voila"),group = GROUPS),
        expand(os.path.join(config['majiq_top_level'],"psi_single",'{sample}' + ".psi.voila"),sample = SAMPLE_NAMES),
        expand(os.path.join(config['majiq_top_level'],"psi_voila_tsv_single",'{sample}' + ".psi.tsv"), sample = SAMPLE_NAMES),
        expand(os.path.join(config['majiq_top_level'],"psi_voila_tsv_single",'{sample}' + "_parsed.csv"), sample = SAMPLE_NAMES),
        os.path.join(config['majiq_top_level'],"psi_voila_tsv_single/" + "full_PSI.csv"),
        expand(os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}" + "_parsed_psi.tsv"),zip, bse = BASES,contrast = CONTRASTS),
        expand(os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}" + "_annotated_junctions.csv"),zip, bse = BASES,contrast = CONTRASTS),
        expand(os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}_annotated.junctions.bed"),zip, bse = BASES,contrast = CONTRASTS)


include: "../rules/majiq_psi.smk"
include: "../rules/majiq_parse_deltas.smk"
include: "../rules/majiq_annotate_junctions.smk"
