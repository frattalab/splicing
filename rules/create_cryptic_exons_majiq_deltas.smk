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

print(SAMPLE_NAMES)


BASES, CONTRASTS = return_bases_and_contrasts()
MAJIQ_DIR = get_output_dir(config['project_top_level'], config['majiq_top_level'])
OUTPUT_CRYPTIC_EXONS = os.path.join(MAJIQ_DIR,"cryptic_exons_beds")
os.system("mkdir -p {0}".format(OUTPUT_CRYPTIC_EXONS))
stringtie_outdir = get_output_dir(config["project_top_level"], config['stringtie_output'])



rule allParse:
    input:
        expand(os.path.join(OUTPUT_CRYPTIC_EXONS,"{bse}_{contrast}_cryptic_exons.bed"),zip, bse = BASES,contrast = CONTRASTS)

rule write_exon_beds:
    input:
        bed = os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv","{bse}_{contrast}_annotated.junctions.bed"),
        assembled_gtf =  os.path.join(stringtie_outdir,"stringtie_merged.gtf")
    output:
        os.path.join(OUTPUT_CRYPTIC_EXONS,"{bse}_{contrast}_cryptic_exons.bed")
    conda:
        "../envs/splicing_dependencies.yml"
    shell:
        """
        Rscript scripts/extract_cryptic_exons_from_gtf.R \
        --transcripts {input.assembled_gtf} \
        --deltabed {input.bed}
        """
