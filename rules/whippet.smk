

import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"
samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))


# had conversation with Ben Whathisface at RNA uk, suggested trying to
# build each bam their own index, and using the PSI for that
# that does mean that you can't use the delte PSI, but let's give her a go

print(SAMPLES)
print(len(SAMPLES))
#our final output
rule all_whippet:
    input:
        expand(config['whippet_output_path'] + "{sample}.sorted.bam.bai", sample = SAMPLES),
		expand(os.path.join(config['whippet_bin'],"index", config['run_name'] + "gencode" + ".jls")),
		expand(os.path.join(config['whippet_bin'],"index", config['run_name'] + "gencode" + ".exons.tab.gz"))




#empty string
bam_dir = config['bam_dir']
bamFiles = ""
for file in os.listdir(bam_dir):
    #for all the aligned sorted bams
    if file.endswith(".Aligned.sorted.out.bam"):
        #{config[samtools_path]} merge takes a list of space separated files
        bamFiles = bamFiles + " " + os.path.join(bam_dir,file)
        quick_temp = os.path.join(bam_dir,file)



rule build_whippet_index:
    input:
        bam = os.path.join(config['bam_dir'], "{sample}" + ".Aligned.sorted.out.bam"),
        bai = os.path.join(config['bam_dir'],  "{sample}" + ".Aligned.sorted.out.bam.bai"),
    output:
        whpt_graph = os.path.join(config['whippet_bin'],"index", "{sample}" + "gencode" + ".jls"),
        whpt_exons = os.path.join(config['whippet_bin'],"index", "{sample}" + "gencode" + ".exons.tab.gz")
    params:
        fasta = config['fasta'],
        gtf = config['gtf'],
        whippet_ind_path = os.path.join(config['whippet_bin'],"whippet-index.jl"),
        annotation_julia_indices = os.path.join(config['whippet_bin'],"index",config['run_name'])
    shell:
        """
        {config[julia]} {params.whippet_ind_path} --fasta {params.fasta} --suppress-low-tsl --bam {input.bam} --gtf {params.gtf} --index {params.annotation_julia_indices}
        """

#first we'll run the psi run
if config['endtype'] == "se":
	rule whippet_psi:
		input:
			config['fastq_files'] + "{sample}.merged.fastq.gz"
		output:
			psi_file = config['whippet_output_path'] + "{sample}.psi.gz",
			jnc_file = config['whippet_output_path'] + "{sample}.jnc.gz",
			sam_file = temp(config['whippet_output_path'] + "{sample}.sam")
		params:
			julia = config['julia'],
			whippet_quant = config['whippet_bin'] + "whippet-quant.jl",
			index = os.path.join(config['whippet_bin'],"index",config['run_name'] + ".jls"),
			output_path = config['whippet_output_path']
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
				psi_file = config['whippet_output_path'] + "{sample}.psi.gz",
				jnc_file = config['whippet_output_path'] + "{sample}.jnc.gz",
				sam_file = temp(config['whippet_output_path'] + "{sample}.sam")
			params:
				julia = config['julia'],
				whippet_quant = config['whippet_bin'] + "whippet-quant.jl",
				index = os.path.join(config['whippet_bin'],"index",config['run_name'] + ".jls"),
				output_path = config['whippet_output_path']
			shell:
				"""
				{params.julia} {params.whippet_quant} {input} --index {params.index} --o {params.output_path} --sam > {output.sam_file}
				"""
else:
	print("End type either 'se' or 'pe'")
#now a quick set of rule to turn all those sam fiels to bam files with samtools; very standard stuff here
rule sam_to_bam:
    input:
        sam_file = config['whippet_output_path'] + "{sample}.sam"
    output:
    #i believe this marks the output file as temporary so that it is deleted one every consuming job has been executed
        temp(config['whippet_output_path'] + "{sample}.bam")
    params:
        samtools = config['samtools_path']
    shell:
        """
        {params.samtools} view -b {input.sam_file} > {output}
        """
rule sort_bams:
    input:
        config['whippet_output_path'] + "{sample}.bam"
    output:
        config['whippet_output_path'] + "{sample}.sorted.bam"
    params:
        samtools = config['samtools_path'],
        old_sam = config['whippet_output_path'] + "{sample}.sam"
    shell:
        """
        rm {params.old_sam}
        {params.samtools} sort {input} > {output}
        """

rule index_sort_bams:
    input:
        config['whippet_output_path'] + "{sample}.sorted.bam"
    output:
        config['whippet_output_path'] + "{sample}.sorted.bam.bai"
    params:
        samtools = config['samtools_path'],
        unsorted_bam = config['whippet_output_path'] + "{sample}.bam"
    shell:
        """
        rm {params.unsorted_bam}
        {params.samtools} index {input} > {output}
        """
