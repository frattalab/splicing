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
#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name']))
BASES, CONTRASTS = return_bases_and_contrasts()
ALLGROUP = list(set(BASES + CONTRASTS))
print(ALLGROUP)


SPECIES = config["species"]
GTF = config['gtf']

#make sure the output folder exists before running anything
bam_dir_temporary = get_output_dir(config["project_top_level"], 'merged_bams')

stringtie_outdir = get_output_dir(config["project_top_level"], config['stringtie_output'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])


#make sure the output folder for STAR exists before running anything
bam_dir = get_output_dir(config["project_top_level"], config['bam_dir'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])

print(scallop_outdir)

rule all_scallop:
    input:
        expand(scallop_outdir + "{grp}.gtf", sample = ALLGROUP),
        os.path.join(scallop_outdir, "scallop_merged.unique.gtf"),
        os.path.join(scallop_outdir, "gffall.scallop_merged.gtf.tmap"),
        expand(os.path.join(scallop_outdir,"{grp}.scallop_merged.gtf"), sample = ALLGROUP)
rule scallop_per_group:
    input:
        bam= os.path.join(bam_dir_temporary,"{grp}.bam"),
        bai= os.path.join(bam_dir_temporary,"{grp}.bam.bai")
    output:
        os.path.join(scallop_outdir,'{grp}' + ".gtf")
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

rule compare_reference:
    input:
        os.path.join(scallop_outdir,'{grp}' + ".gtf")
    output:
        os.path.join(scallop_outdir, '{grp}' + "gffall.scallop_merged.gtf.tmap")
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
