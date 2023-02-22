import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "../rules/helpers.py"

configfile: "config/config.yaml"

#reading in the samples and dropping the samples to be excluded in order to get a list of sample names

BASES, CONTRASTS = return_bases_and_contrasts()
MAJIQ_DIR = get_output_dir(config['project_top_level'], config['majiq_top_level'])

samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name']))
BASES, CONTRASTS = return_bases_and_contrasts()
ALLGROUP = list(set(BASES + CONTRASTS))
print(ALLGROUP)


SPECIES = config["species"]
GTF = config['gtf']

#make sure the output folder for STAR exists before running anything
bam_dir = get_output_dir(config["project_top_level"], config['bam_dir'])
bam_dir_temporary = get_output_dir(config["project_top_level"], 'merged_bams')

stringtie_outdir = get_output_dir(config["project_top_level"], config['stringtie_output'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])

both_output_dirs = [stringtie_outdir,scallop_outdir]
def return_awk():
    awk_string = """awk '{if ($7 != ".") {print}}'"""
    return awk_string

awk_string = """awk '{if ($7 != ".") {print}}'"""

print(both_output_dirs)
print(awk_string)
rule allMerging:
    input:
        expand(os.path.join(scallop_outdir,'{grp}' + ".gtf"),grp = ALLGROUP),
        expand(os.path.join(stringtie_outdir,'{grp}' + ".gtf"),grp = ALLGROUP),
        expand(os.path.join(stringtie_outdir,"{grp}.unstranded_filtered.unique.gtf"),grp = ALLGROUP),
        expand('{outputdir}{grp}' + ".annotated.gtf",outputdir =both_output_dirs, grp = ALLGROUP),
        expand("{outputdir}" + "{grp}.unique.gtf",outputdir =both_output_dirs, grp = ALLGROUP),
        return_bed_and_bases(BASES,CONTRASTS,both_output_dirs)


rule merge_bam_groups:
    input:
        group_bam_files = lambda wildcards: file_path_list(wildcards.grp,bam_dir,config['bam_suffix'] + '.bam')
    output:
        bam= temp(os.path.join(bam_dir_temporary,"{grp}.bam"))
    wildcard_constraints:
        grp="|".join(ALLGROUP)
    threads: 10
    shell:
        """
        samtools merge {output.bam} {input.group_bam_files} --threads {threads}
        """

rule index_merged_group:
    input:
        bam = os.path.join(bam_dir_temporary,"{grp}.bam")
    output:
        bai= temp(os.path.join(bam_dir_temporary,"{grp}.bam.bai"))
    wildcard_constraints:
        grp="|".join(ALLGROUP)
    threads: 10
    shell:
        """
        samtools index {input.bam} {output.bai}
        """

rule scallop_per_group:
    input:
        bam= os.path.join(bam_dir_temporary,"{grp}.bam"),
        bai= os.path.join(bam_dir_temporary,"{grp}.bam.bai")
    output:
        os.path.join(scallop_outdir,'{grp}' + ".gtf")
    wildcard_constraints:
        grp="|".join(ALLGROUP)
    conda:
        "../envs/scallop.yml"
    params:
        scallop_path = 'scallop2',
        scallop_out_folder = scallop_outdir,
        scallop_extra_config = return_parsed_extra_params(config['scallop_extra_parameters']),
        libtype = config['scallop_strand']
    shell:
        """
        mkdir -p {params.scallop_out_folder}
        {params.scallop_path} \
        -i {input.bam} \
        -o {output} \
        --library_type {params.libtype} \
        {params.scallop_extra_config}
        """

rule stringtie_per_group:
    input:
        bam= os.path.join(bam_dir_temporary,"{grp}.bam"),
        bai= os.path.join(bam_dir_temporary,"{grp}.bam.bai")
    output:
        stringtie_outdir + "{grp}.gtf"
    wildcard_constraints:
        grp="|".join(ALLGROUP)
    conda:
        "../envs/stringtie.yml"
    shell:
        "stringtie {input.bam} -G {GTF} -o {output}"


rule compare_reference:
    input:
        "{outputdir}" + "{grp}.gtf"
    output:
        "{outputdir}{grp}.annotated.gtf",
        "{outputdir}" + "{grp}.{grp}.gtf.tmap"
    wildcard_constraints:
        grp="|".join(ALLGROUP)
    params:
        ref_gtf = GTF,
        gffcompare = config['gffcompare'],
        prefix = os.path.join('{outputdir}', '{grp}')
    shell:
        """
        {params.gffcompare} -o {params.prefix} -r {params.ref_gtf} {input}
        """

rule fetch_unique:
    input:
        sample_gtf = "{outputdir}" + "{grp}.gtf",
        sample_tmap = "{outputdir}" + "{grp}.{grp}.gtf.tmap"
    output:
        "{outputdir}" + "{grp}.unique.gtf"
    wildcard_constraints:
        grp="|".join(ALLGROUP),
        outputdir="|".join(both_output_dirs)
    params:
        ref_gtf = GTF,
        gtfcuff = config['gtfcuff']
    shell:
        """
        {params.gtfcuff} puniq {input.sample_tmap} {input.sample_gtf} {params.ref_gtf} {output}
        """

rule filter_stringtie:
    #stringtie outputs strandless 1 exon transcripts, R crashes when it reads those in, so filter them out
    input:
        stringtie_outdir + "{grp}.unique.gtf"
    output:
        stringtie_outdir + "{grp}.unstranded_filtered.unique.gtf"
    params:
        my_awk = lambda wildcards: return_awk()
    shell:
        """
        {params.my_awk} {input} > {output}
        """
rule filter_scallop:
    #stringtie outputs strandless 1 exon transcripts, R crashes when it reads those in, so filter them out
    input:
        scallop_outdir + "{grp}.unique.gtf"
    output:
        scallop_outdir + "{grp}.unstranded_filtered.unique.gtf"
    shell:
        """
        ln -s {input} {output}
        """
rule write_exon_beds:
    input:
        delta_csv = os.path.join(MAJIQ_DIR,"delta_psi_voila_tsv","{bse}-{contrast}_annotated_junctions.csv"),
        assembled_gtf =  "{outputdir}" + "{contrast}.unstranded_filtered.unique.gtf"
        # assembled_gtf =  os.path.join(stringtie_outdir,"stringtie_merged.gtf")
    output:
        "{outputdir}" + "{bse}-{contrast}_cryptic_exons.bed",
    # conda:
    #     "../envs/splicing_dependencies.yml"
    wildcard_constraints:
        outputdir="|".join(both_output_dirs)
    shell:
        """
        Rscript scripts/extract_cryptic_exons_from_gtf.R \
        --transcripts {input.assembled_gtf} \
        --delta {input.delta_csv} \
        --outputname {output}
        """