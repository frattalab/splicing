import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "helpers.py"
localrules: compose_gtf_list
##############################
##### STOLEN FROM https://github.com/bioinformatics-core-shared-training/RNAseq_March_2019/tree/master/
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
stringtie_outdir = get_output_dir(config["project_top_level"], config['stringtie_output'])

print(stringtie_outdir)

rule all_stringtie:
    input:
        expand(stringtie_outdir + "{sample}.assemble.gtf", sample = SAMPLE_NAMES),
        expand(stringtie_outdir + "gffall.{sample}.gtf.tmap",sample = SAMPLE_NAMES),
        expand(stringtie_outdir + "{sample}.unique.gtf",sample = SAMPLE_NAMES),
        os.path.join(stringtie_outdir,"scallop_merged.gtf")

rule StringTie_Assemble:
    input:
        bam = lambda wildcards: star_outdir + '{sample}' + config['bam_suffix'],
        ref_gtf = GTF
    output:
        stringtie_outdir + "{sample}.assemble.gtf"
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie {input.bam} -G {input.ref_gtf} -o {output}"


rule compare_reference:
    input:
        os.path.join(stringtie_outdir + "{sample}.assemble.gtf")
    output:
        os.path.join(stringtie_outdir, "gffall.{sample}.gtf.tmap")
    params:
        ref_gtf = GTF,
        gffcompare = config['gffcompare']
    shell:
        """
        {params.gffcompare} -o gffall -r {params.ref_gtf} {input}
        """

rule fetch_unique:
    input:
        sample_tmap = os.path.join(stringtie_outdir, "gffall.{sample}.gtf.tmap"),
        sample_gtf = os.path.join(stringtie_outdir,'{sample}' + ".assemble.gtf")
    output:
        os.path.join(stringtie_outdir, "{sample}.unique.gtf")
    params:
        ref_gtf = GTF,
        gtfcuff = config['gtfcuff']
    shell:
        """
        {params.gtfcuff} puniq {input.sample_tmap} {input.sample_gtf} {params.ref_gtf} {output}
        """

rule compose_gtf_list:
    input:
        expand(stringtie_outdir + "{sample}.unique.gtf", sample=SAMPLE_NAMES)
    output:
        txt = os.path.join(stringtie_outdir,"gtf_list.txt")
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)

rule merge_stringtie_gtfs:
    input:
        gtf_list = os.path.join(stringtie_outdir,"gtf_list.txt")
    output:
        merged_gtf = os.path.join(stringtie_outdir,"scallop_merged.gtf")
    params:
        gtfmerge = config['gtfmerge']
    shell:
        """
        {params.gtfmerge} union {input.gtf_list} {output.merged_gtf} -t 2 -n
        """
