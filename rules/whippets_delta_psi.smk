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

whippet_folder = os.path.join(config['top_level_project_folder'],config['whippet_output_path'])
whippet_psi_output_folder = os.path.join(config['top_level_project_folder'],config['whippet_output_path'],'psi/')
rule allWhippetDeltaPSI:
    input:
        expand(whippet_folder + "deltapsi/","{bse}_{contrast}" + ".diff.gz"),zip, bse = BASES,contrast = CONTRASTS)

rule whippet_delta_psi:
    input:
        base_group_whipp = lambda wildcards: whippets_psi_files_from_contrast(wildcards.bse,whippet_psi_output_folder),
        contrast_group_whipp = lambda wildcards: whippets_psi_files_from_contrast(wildcards.contrast,whippet_psi_output_folder)
    output:
        whippet_folder + "deltapsi/" + "{bse}_{contrast}" + ".diff.gz"
    params:
        julia = config['julia'],
        whippet_quant = config['whippet_bin'] + "whippet-delta.jl",
        output_path = whippet_folder + "deltapsi/" + "{bse}_{contrast}"
    shell:
        """
        export JULIA_PKGDIR=/SAN/vyplab/alb_projects/tools/julia_pkgdir/v0.6/
        mkdir -p {params.delta_psi_output_folder}
        {params.julia} {params.whippet_quant}
        -a {input.base_group_whipp}, -b {input.contrast_group_whipp}, \
        --out {params.output_path}
        """
