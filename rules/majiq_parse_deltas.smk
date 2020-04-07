import pandas as pd
import os
import subprocess
import yaml

configfile: "config/config.yaml"
include: "helpers.py"
include: "scripts/helpers.py"
localrules: create_majiq_config_file
#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))

BASES, CONTRASTS = return_bases_and_contrasts()
print(BASES)
print(CONTRASTS)

rule all:
    input:
        # expand(config['majiq_top_level'] + "delta_psi/" + "{bse}_{contrast}" + ".deltapsi.tsv",zip, bse = BASES,contrast = CONTRASTS),
        # expand(os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}" + ".psi.tsv"),zip, bse = BASES,contrast = CONTRASTS),
        # config['majiq_top_level'] + config['run_name'] + "_majiqConfig.tsv",
        # expand(os.path.join(config['majiq_top_level'],"builder",'{name}' + ".majiq"),name = SAMPLE_NAMES),
        # expand(os.path.join(config['majiq_top_level'],"psi",'{group}' + ".psi.voila"),group = GROUPS),
        expand(os.path.join(config['majiq_top_level'],"star_beds",'{name}' + ".bed"),name = SAMPLE_NAMES)

rule star_tabs_to_beds:
    input:
        sj_tab = lambda wildcards: os.path.join(config['bam_dir'],"{name}" + ".SJ.out.tab")
    output:
        bed = os.path.join(config['majiq_top_level'],"star_beds",'{name}' + ".bed")
    shell:
        """
        mkdir -p {config['majiq_top_level']}/star_beds/
        python3 ./scripts/splicejunction2bed.py -i {input.sj_tab} -o {output.bed}
        """
# rule majiq_delta_psi_tsv:
#     input:
#     #this is always calling from the column named 'group' in the sample csv file
#         voila_file = lambda wildcards: os.path.join(config['majiq_top_level'],"delta_psi","{bse}_{contrast}" + ".deltapsi.voila")
#     output:
#         tsv = os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}" + ".psi.tsv")
#     params:
#         voila_path = config['voila_path'],
#         psi_output_folder = os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv"),
#         splice_graph = os.path.join(config['majiq_top_level'],"builder", "splicegraph.sql"),
#         extra_voila_paramters = return_parsed_extra_params(config['extra_voila_parameters'])
#     shell:
#         """
#         mkdir -p {params.psi_output_folder}
#         {params.voila_path} tsv {params.splice_graph} {input.voila_file} -f {output.tsv} {params.extra_voila_paramters}
#         """
