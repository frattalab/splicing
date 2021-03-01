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
        expand(os.path.join(MAJIQ_DIR,"psi_voila_tsv_single",'{sample}' + "_parsed.csv"), sample = SAMPLE_NAMES),
        expand(os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv","{bse}_{contrast}" + "_parsed_psi.tsv"),zip, bse = BASES,contrast = CONTRASTS),
        # os.path.join(MAJIQ_DIR,"psi_voila_tsv_single/" + "full_PSI.csv"),

rule write_exon_beds:
    input:
        bed = os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv","{bse}_{contrast}_annotated.junctions.bed"),
        assembled_gtf =  os.path.join(config["project_top_level"],"all_assemblers_merged.gtf")
    output:
        os.path.join(MAJIQ_DIR,,"{bse}_{contrast}_cryptic_exons.bed")
    conda:
        "../envs/splicing_dependencies.yml"
    params:
        extra_junction_parameters = return_parsed_extra_params(config['annotated_junctions_extra_parameters']),
        trackname = "{bse}_{contrast}.bed"
    shell:
        """
        Rscript scripts/make_bed_from_annotated_command_line.R \
        --parsed {input.csv} \
        --out {output} \
        {params.extra_junction_parameters}
        """
