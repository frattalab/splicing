import pandas as pd
import os
import subprocess
import yaml

configfile: "config/config.yaml"

include: "helpers.py"

#reading in the samples and dropping the samples to be excluded in order to get a list of sample names

samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
SJ_NAMES  = list(set(samples2['sample_name']))

GROUPS = list(set(samples2['group']))

BASES, CONTRASTS = return_bases_and_contrasts()

rule all:
    input:
        expand(os.path.join(config['majiq_top_level'],"psi_voila_tsv_single",'{sample}' + "_parsed.csv"), sample = SAMPLE_NAMES),
        os.path.join(config['majiq_top_level'],"psi_voila_tsv_single/" + "full_PSI.csv")
        # expand(os.path.join(config['majiq_top_level'],"star_beds",'{sjname}' + ".bed"),sjname = SJ_NAMES)

# rule star_tabs_to_beds:
#     input:
#         sj_tab = os.path.join(config['bam_dir'],"{sjname}" + ".SJ.out.tab")
#     params:
#         out_folder = config['majiq_top_level'] + "star_beds/"
#     output:
#         bed = os.path.join(config['majiq_top_level'],"star_beds",'{sjname}' + ".bed")
#     shell:
#         """
#         mkdir -p {params.out_folder}
#         python3 ./scripts/splicejunction2bed.py -i {input.sj_tab} -o {output.bed}
#         """
rule majiq_psi_parse:
    input:
        psi_voila_tsv = lambda wildcards: os.path.join(config['majiq_top_level'],"psi_voila_tsv_single",'{sample}' + ".psi.tsv")
    wildcard_constraints:
        name="|".join(SAMPLE_NAMES)
    output:
        parsed_csv = os.path.join(config['majiq_top_level'],"psi_voila_tsv_single",'{sample}' + "_parsed.csv")
    shell:
        """
        Rscript scripts/parsing_psi_command_line.R --input {input.psi_voila_tsv} -o {output.parsed_csv}
        """

rule combine_psi_per_sample:
    input:
        all_parsed_csvs = get_single_psi_parsed_files()
    output:
        parsed_csv = os.path.join(config['majiq_top_level'],"psi_voila_tsv_single/" + "full_PSI.csv")
    params:
        psi_folder = os.path.join(config['majiq_top_level'],"psi_voila_tsv_single")
    shell:
        """
        Rscript scripts/writing_final_parsed_psi_command_line.R --folder {params.psi_folder} --out {output.parsed_csv}
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
