

import os,re

configfile: "config/config.yaml"
#lazy way of pulling out all the sample names in the fastq fiels
if config['endtype'] == "se":
	SAMPLES = [re.sub(".merged.fastq.gz", "",smp) for smp in os.listdir(config['fastq_files'])]
elif config['endtype'] == "pe":
	SAMPLES = [re.sub("_\d.merged.fastq.gz", "",smp) for smp in os.listdir(config['fastq_files'])]
else:
	print("End type either 'se' or 'pe'")
import pandas as pd
import os
import subprocess

#first rule of a snakemake file always defines the final desired output. Eventually this file should produced the psi file
#and the sorted bam
print(SAMPLES)
print(len(SAMPLES))
#our final output
rule all_whippet:
    input:
        expand(config['data_output_path'] + "{sample}.sorted.bam.bai", sample = SAMPLES)

#first we'll run the psi run
if config['endtype'] == "se":
	rule whippet_psi:
		input:
			config['fastq_files'] + "{sample}.merged.fastq.gz"
		output:
			psi_file = config['data_output_path'] + "{sample}.psi.gz",
			jnc_file = config['data_output_path'] + "{sample}.jnc.gz",
			sam_file = temp(config['data_output_path'] + "{sample}.sam")
		params:
			julia = config['julia'],
			whippet_quant = config['whippet_bin'] + "whippet-quant.jl",
			index = os.path.join(config['whippet_bin'],"index",config['run_name'] + ".jls"),
			output_path = config['data_output_path']
		shell:
			"""
			{params.julia} {params.whippet_quant} {input} --index {params.index} --o {params.output_path} --sam > {output.sam_file}
			"""
elif config['endtype'] == "pe":
		rule whippet_psi:
			input:
				fwd = config['fastq_files'] + "{sample}_1.merged.fastq.gz",
				rev = config['fastq_files'] + "{sample}_2.merged.fastq.gz"
			output:
				psi_file = config['data_output_path'] + "{sample}.psi.gz",
				jnc_file = config['data_output_path'] + "{sample}.jnc.gz",
				sam_file = temp(config['data_output_path'] + "{sample}.sam")
			params:
				julia = config['julia'],
				whippet_quant = config['whippet_bin'] + "whippet-quant.jl",
				index = os.path.join(config['whippet_bin'],"index",config['run_name'] + ".jls"),
				output_path = config['data_output_path']
			shell:
				"""
				{params.julia} {params.whippet_quant} {input} --index {params.index} --o {params.output_path} --sam > {output.sam_file}
				"""
else:
	print("End type either 'se' or 'pe'")
#now a quick set of rule to turn all those sam fiels to bam files with samtools; very standard stuff here
rule sam_to_bam:
    input:
        sam_file = config['data_output_path'] + "{sample}.sam"
    output:
    #i believe this marks the output file as temporary so that it is deleted one every consuming job has been executed
        temp(config['data_output_path'] + "{sample}.bam")
    params:
        samtools = config['samtools_path']
    shell:
        """
        {params.samtools} view -b {input.sam_file} > {output}
        """
rule sort_bams:
    input:
        config['data_output_path'] + "{sample}.bam"
    output:
        config['data_output_path'] + "{sample}.sorted.bam"
    params:
        samtools = config['samtools_path'],
        old_sam = config['data_output_path'] + "{sample}.sam"
    shell:
        """
        rm {params.old_sam}
        {params.samtools} sort {input} > {output}
        """

rule index_sort_bams:
    input:
        config['data_output_path'] + "{sample}.sorted.bam"
    output:
        config['data_output_path'] + "{sample}.sorted.bam.bai"
    params:
        samtools = config['samtools_path'],
        unsorted_bam = config['data_output_path'] + "{sample}.bam"
    shell:
        """
        rm {params.unsorted_bam}
        {params.samtools} index {input} > {output}
        """
