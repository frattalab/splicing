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

rule allAnnotated:
    input:
        expand(os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}" + "_annotated_junctions.gff3"),zip, bse = BASES,contrast = CONTRASTS),
        expand(os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}_annotated.junctions.bed"),zip, bse = BASES,contrast = CONTRASTS)

rule annotatate_delta:
    input:
        tsv = os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}" + "_parsed_psi.tsv")
    output:
        os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}_annotated_junctions.gff3")
    params:
        gtf = config['gtf'],
        psi_output_folder = os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}_annotated_junctions")
    shell:
        """
        Rscript scripts/add_junction_annotations_command_line.R --deltapsi {input.tsv} --out {params.psi_output_folder} --gtf {params.gtf}
        """

rule write_junctions_beds:
    input:
        csv = os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}_annotated_junctions.csv")
    output:
        os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}_annotated.junctions.bed")
    params:
        extra_junction_parameters = return_parsed_extra_params(config['annotated_junctions_extra_parameters']),
        psi_output_folder = os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}_annotated_junctions"),
        trackname = "{bse}_{contrast}.bed"
    shell:
        """
        Rscript scripts/make_bed_from_annotated_command_line.R --parsed {input.csv} \
        --out {params.psi_output_folder} \
        {params.extra_junction_parameters}
        """
