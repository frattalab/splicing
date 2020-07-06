import pandas as pd
import os
import subprocess
import yaml
# note to self - probaly should not be hard coding threads like i do
configfile: "config/config.yaml"
include: "helpers.py"

#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))

BASES, CONTRASTS = return_bases_and_contrasts()
print(BASES)
print(CONTRASTS)

whippet_delta_psi_output_folder = os.path.join(config['top_level_project_folder'],\
                                               config['whippet_output_path'],"delta_psi/")

rule allPSI:
    input:
        expand(whippet_delta_psi_output_folder,"{bse}_{contrast}" + ".diff.gz"),zip, bse = BASES,contrast = CONTRASTS)

rule whippet_delta_psi:
    input:
        base_group_majiq = lambda wildcards: whippets_psi_files_from_contrast(wildcards.bse,top_level_folder),
        contrast_group_majiq = lambda wildcards: whippets_psi_files_from_contrast(wildcards.contrast,top_level_folder)
    output:
        os.path.join(config['majiq_top_level'],"delta_psi","{bse}_{contrast}" + ".deltapsi.tsv"),
        os.path.join(config['majiq_top_level'],"delta_psi","{bse}_{contrast}" + ".deltapsi.voila")
    params:
        majiq_path = config['majiq_path'],
        delta_psi_output_folder = os.path.join(config['majiq_top_level'],"delta_psi"),
        majiq_psi_extra_parameters = return_parsed_extra_params(config['majiq_psi_extra_parameters'])
    threads:
        8
    shell:
        """
        mkdir -p {params.delta_psi_output_folder}
        {params.majiq_path} deltapsi -grp1 {input.base_group_majiq} -grp2 {input.contrast_group_majiq} -j {threads} -o {params.delta_psi_output_folder}{params.majiq_psi_extra_parameters} --name {wildcards.bse} {wildcards.contrast}
        """
