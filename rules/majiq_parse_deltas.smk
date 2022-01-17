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
SJ_NAMES  = list(set(samples2['sample_name']))
print(SAMPLE_NAMES)
GROUPS = list(set(samples2['group']))

BASES, CONTRASTS = return_bases_and_contrasts()
MAJIQ_DIR = get_output_dir(config['project_top_level'], config['majiq_top_level'])

rule allParse:
    input:
        # expand(os.path.join(MAJIQ_DIR,"psi_voila_tsv_single",'{sample}' + "_parsed.csv"), sample = SAMPLE_NAMES),
        expand(os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv","{bse}-{contrast}" + "_parsed_psi.tsv"),zip, bse = BASES,contrast = CONTRASTS),
        # os.path.join(MAJIQ_DIR,"psi_voila_tsv_single/" + "full_PSI.csv"),


        # expand(os.path.join(MAJIQ_DIR,"star_beds",'{sjname}' + ".bed"),sjname = SJ_NAMES)

# rule star_tabs_to_beds:
#     input:
#         sj_tab = os.path.join(config['bam_dir'],"{sjname}" + ".SJ.out.tab")
#     params:
#         out_folder = MAJIQ_DIR + "star_beds/"
#     output:
#         bed = os.path.join(MAJIQ_DIR,"star_beds",'{sjname}' + ".bed")
#     shell:
#         """
#         mkdir -p {params.out_folder}
#         python3 ./scripts/splicejunction2bed.py -i {input.sj_tab} -o {output.bed}
#         """
rule majiq_psi_parse:
    input:
        psi_voila_tsv = lambda wildcards: os.path.join(MAJIQ_DIR,"psi_voila_tsv_single",'{sample}' + ".psi.tsv")
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    conda:
        "../envs/splicing_dependencies.yml"
    output:
        parsed_csv = os.path.join(MAJIQ_DIR,"psi_voila_tsv_single",'{sample}' + "_parsed.csv")
    shell:
        """
        Rscript scripts/parsing_psi_command_line.R --input {input.psi_voila_tsv} -o {output.parsed_csv}
        """

# rule combine_psi_per_sample:
#     input:
#         all_parsed_csvs = get_single_psi_parsed_files()
#     output:
#         parsed_csv = os.path.join(MAJIQ_DIR,"psi_voila_tsv_single/" + "full_PSI.csv")
#     params:
#         psi_folder = os.path.join(MAJIQ_DIR,"psi_voila_tsv_single")
#     shell:
#         """
#         Rscript scripts/writing_final_parsed_psi_command_line.R --folder {params.psi_folder} --out {output.parsed_csv}
#         """
#

rule majiq_finish_deltas:
    input:
        tsv = lambda wildcards: os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv","{bse}-{contrast}" + ".psi.tsv")
    output:
        parsed_tsv = os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv","{bse}-{contrast}" + "_parsed_psi.tsv")
    conda:
        "../envs/splicing_dependencies.yml"
    params:
        psi_output_folder = os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv")
    shell:
        """
        mkdir -p {params.psi_output_folder}
        Rscript scripts/write_voila_delta_parsed.R --deltafile {input.tsv} --out {output.parsed_tsv}
        """
