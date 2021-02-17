import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "helpers.py"
localrules: compose_gtf_list
##############################
##### Final output is a merged gtf of all samples containing transcripts that
##### were found in the samples but NOT the input GTF
##############################
#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name']))

print(SAMPLE_NAMES)

SPECIES = config["species"]
GTF = get_gtf(SPECIES)

#make sure the output folder for STAR exists before running anything
star_outdir = get_output_dir(config["project_top_level"], config['star_output_folder'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])
print(scallop_outdir)
rule all_scallop:
    input:
        expand(scallop_outdir + '{sample}' + ".gtf", sample = SAMPLE_NAMES),
        expand(scallop_outdir + "gffall.{sample}.gtf.tmap",sample = SAMPLE_NAMES),
        expand(scallop_outdir + "{sample}.unique.gtf",sample = SAMPLE_NAMES),
        os.path.join(scallop_outdir,"scallop_merged.gtf")


rule scallop_per_samp:
    input:
        bam_file = lambda wildcards: star_outdir + '{sample}' + config['bam_suffix']
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    output:
        os.path.join(scallop_outdir,'{sample}' + ".gtf")
    params:
        scallop_path = config['scallop_path'],
        verbose = 0,
        scallop_out_folder = scallop_outdir,
        scallop_extra_config = return_parsed_extra_params(config['scallop_extra_parameters']),
        libtype = get_scallop_strand(config["feature_counts_strand_info"]),
    shell:
        """
        mkdir -p {params.scallop_out_folder}
        {params.scallop_path} \
        -i {input.bam_file} \
        -o {output} \
        --library_type {params.libtype} \
        --verbose {params.verbose} \
        {params.scallop_extra_config}
        """

rule compare_reference:
    input:
        os.path.join(scallop_outdir,'{sample}' + ".gtf")
    output:
        os.path.join(scallop_outdir, "gffall.{sample}.gtf.tmap")
    params:
        ref_gtf = GTF,
        gffcompare = config['gffcompare']
    shell:
        """
        {params.gffcompare} -o gffall -r {params.ref_gtf} {input}
        """

rule fetch_unique:
    input:
        sample_tmap = os.path.join(scallop_outdir, "gffall.{sample}.gtf.tmap"),
        sample_gtf = os.path.join(scallop_outdir,'{sample}' + ".gtf")
    output:
        os.path.join(scallop_outdir, "{sample}.unique.gtf")
    params:
        ref_gtf = GTF,
        gtfcuff = config['gtfcuff']
    shell:
        """
        {params.gtfcuff} puniq {input.sample_tmap} {input.sample_gtf} {params.ref_gtf} {output}
        """

rule compose_gtf_list:
    input:
        expand(scallop_outdir + "{sample}.unique.gtf", sample=SAMPLE_NAMES)
    output:
        txt = os.path.join(scallop_outdir,"gtf_list.txt")
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)

rule merge_scallop_gtfs:
    input:
        gtf_list = os.path.join(scallop_outdir,"gtf_list.txt")
    output:
        merged_gtf = os.path.join(scallop_outdir,"scallop_merged.gtf")
    params:
        gtfmerge = config['gtfmerge']
    shell:
        """
        {params.gtfmerge} union {input.gtf_list} {output.merged_gtf} -t 2 -n
        """
