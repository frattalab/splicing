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

BASES, CONTRASTS = return_bases_and_contrasts()
print(BASES)
print(CONTRASTS)

rule allAnnotated:
    input:
        expand(os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}" + "_annotated_junctions.gff3"),zip, bse = BASES,contrast = CONTRASTS)


rule annotatate_delta:
    input:
        os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}_parsed_psi.tsv"
    output:
        os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}_annotated_junctions.gff3")
    params:
        gtf = config['gtf']
    shell:
        """
        Rscript scripts/add_junction_annotations_command_line.R --deltapsi {input.tsv} --out {params.psi_output_folder} --gtf {params.gtf}
        """
