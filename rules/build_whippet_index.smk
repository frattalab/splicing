import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"



#empty string
bam_dir = config['bam_dir']
bamFiles = ""
for file in os.listdir(bam_dir):
    #for all the aligned sorted bams
    if file.endswith(".Aligned.sorted.out.bam"):
        #{config[samtools_path]} merge takes a list of space separated files
        bamFiles = bamFiles + " " + os.path.join(bam_dir,file)
        quick_temp = os.path.join(bam_dir,file)


#final output of the index run is a graph and an exons tab
rule all_index:
    input:
        os.path.join(config['whippet_bin'],"index", config['run_name'] + "gencode" + ".jls"),
        os.path.join(config['whippet_bin'],"index", config['run_name'] + "gencode" + ".exons.tab.gz")

#just a fast rule to merge the bams produced previously
rule merge_bam:
    input:
        quick_temp
    output:
        output_bam = temp(os.path.join(config['data_output_path'],config['run_name'] + "_merged.bam"))
    shell:
        """
        {config[samtools_path]} merge {output.output_bam}{bamFiles}
        """
#next three rules are pretty simple, sort, remove dups, index with samtools
rule sort_bam:
    input:
        merged_bam = os.path.join(config['data_output_path'],config['run_name'] + "_merged.bam")
    output:
        sorted_mergedbam = temp(os.path.join(config['data_output_path'], config['run_name'] + "_merged.sorted.bam"))
    shell:
        """
        {config[samtools_path]} sort -o {output.sorted_mergedbam} {input.merged_bam}
        """
rule remove_dups:
    input:
        sorted_mergedbam = os.path.join(config['data_output_path'], config['run_name'] + "_merged.sorted.bam")
    output:
        sorted_rmdupsbam = temp(os.path.join(config['data_output_path'], config['run_name'] + "_merged.sorted.rmdup.bam"))
    shell:
        """
        {config[samtools_path]} rmdup {input.sorted_mergedbam} {output.sorted_rmdupsbam}
        """
rule index_bam:
    input:
        sorted_rmdupsbam = os.path.join(config['data_output_path'], config['run_name'] + "_merged.sorted.rmdup.bam")
    output:
        sorted_rmdupsbambai = temp(os.path.join(config['data_output_path'], config['run_name'] + "_merged.sorted.rmdup.bam.bai"))
    shell:
        """
        {config[samtools_path]} index {input.sorted_rmdupsbam}
        """

rule build_whippet_index:
    input:
        bam = os.path.join(config['data_output_path'], config['run_name'] + "_merged.sorted.rmdup.bam"),
        bai = os.path.join(config['data_output_path'], config['run_name'] + "_merged.sorted.rmdup.bam.bai"),
    output:
        whpt_graph = os.path.join(config['whippet_bin'],"index", config['run_name'] + "gencode" + ".jls"),
        whpt_exons = os.path.join(config['whippet_bin'],"index", config['run_name'] + "gencode" + ".exons.tab.gz")
    params:
        fasta = config['fasta'],
        gtf = config['gtf'],
        whippet_ind_path = os.path.join(config['whippet_bin'],"whippet-index.jl"),
        annotation_julia_indices = os.path.join(config['whippet_bin'],"index",config['run_name'])
    shell:
        """
        {config[julia]} {params.whippet_ind_path} --fasta {params.fasta} --suppress-low-tsl --bam {input.bam} --gtf {params.gtf} --index {params.annotation_julia_indices}
        """
