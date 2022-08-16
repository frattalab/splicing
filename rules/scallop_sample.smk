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
BASES, CONTRASTS = return_bases_and_contrasts()


GTF = config['gtf']

#make sure the output folder for STAR exists before running anything
bam_dir = get_output_dir(config["project_top_level"], config['bam_dir'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])

print(scallop_outdir)

rule all_scallop:
    input:
        expand(scallop_outdir + "{sample}.gtf", sample = SAMPLE_NAMES),
        os.path.join(scallop_outdir, "scallop_merged.unique.gtf"),
        os.path.join(scallop_outdir, "gffall.scallop_merged.gtf.tmap"),
        expand(os.path.join(scallop_outdir,"{bse}.scallop_merged.gtf"), bse = BASES),
        expand(os.path.join(scallop_outdir,"{contrast}.scallop_merged.gtf"), contrast = CONTRASTS)

rule scallop_per_samp:
    input:
        bam_file = lambda wildcards: bam_dir + '{sample}' + config['bam_suffix'] + ".bam"
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    output:
        os.path.join(scallop_outdir,'{sample}' + ".gtf")
    conda:
        "../envs/scallop.yml"
    params:
        scallop_path = 'scallop2',
        verbose = 0,
        scallop_out_folder = scallop_outdir,
        scallop_extra_config = return_parsed_extra_params(config['scallop_extra_parameters']),
        libtype = config['scallop_strand']
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

rule compose_gtf_list_everyone:
    input:
        expand(scallop_outdir + "{sample}.gtf", sample=SAMPLE_NAMES)
    output:
        txt = os.path.join(scallop_outdir,"gtf_list.txt")
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)

rule merge_scallop_gtfs_everyone:
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
######### DOING THE SAME THING THREE TIMES NOW
rule compose_gtf_list_bases:
    input:
        base_group_scallop = lambda wildcards: file_path_list(wildcards.bse,scallop_outdir,".gtf")
    wildcard_constraints:
        bse="|".join(BASES)
    output:
        txt = temp(os.path.join(scallop_outdir,"{bse}.gtf_list.txt"))
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)

rule merge_scallop_gtfs_bases:
    input:
        gtf_list = os.path.join(scallop_outdir,"{bse}.gtf_list.txt")
    output:
        merged_gtf = os.path.join(scallop_outdir,"{bse}.scallop_merged.gtf")
    wildcard_constraints:
        bse="|".join(BASES)
    params:
        gtfmerge = config['gtfmerge']
    shell:
        """
        {params.gtfmerge} union {input.gtf_list} {output.merged_gtf} -t 2 -n
        """

rule compose_gtf_list_contrast:
    input:
        base_group_scallop = lambda wildcards: file_path_list(wildcards.contrast,scallop_outdir,".gtf")
    wildcard_constraints:
        contrast="|".join(CONTRASTS)
    output:
        txt = temp(os.path.join(scallop_outdir,"{contrast}.gtf_list.txt"))
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)

rule merge_scallop_gtfs_contrasts:
    input:
        gtf_list = os.path.join(scallop_outdir,"{contrast}.gtf_list.txt")
    wildcard_constraints:
        contrast="|".join(CONTRASTS)
    output:
        merged_gtf = os.path.join(scallop_outdir,"{contrast}.scallop_merged.gtf")
    params:
        gtfmerge = config['gtfmerge']
    shell:
        """
        {params.gtfmerge} union {input.gtf_list} {output.merged_gtf} -t 2 -n
        """

rule compare_reference:
    input:
        os.path.join(scallop_outdir,"scallop_merged.gtf")
    output:
        os.path.join(scallop_outdir, "gffall.scallop_merged.gtf.tmap")
    params:
        ref_gtf = GTF,
        gffcompare = config['gffcompare']
    shell:
        """
        {params.gffcompare} -o gffall -r {params.ref_gtf} {input}
        """

rule fetch_unique:
    input:
        sample_tmap = os.path.join(scallop_outdir,"scallop_merged.gtf"),
        sample_gtf = os.path.join(scallop_outdir, "gffall.scallop_merged.gtf.tmap")
    output:
        os.path.join(scallop_outdir, "scallop_merged.unique.gtf")
    params:
        ref_gtf = GTF,
        gtfcuff = config['gtfcuff']
    shell:
        """
        {params.gtfcuff} puniq {input.sample_tmap} {input.sample_gtf} {params.ref_gtf} {output}
        """
