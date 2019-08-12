import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"


#empty string
bams = ""
bam_dir = config['bam_dir']
for file in os.listdir(bam_dir):
	#for all the aligned sorted bams
    if file.endswith(".Aligned.sorted.out.bam"):
    	#samtools merge takes a list of space separated files
        bams = bams + " " + os.path.join(bam_dir,file)

#final output of the index run is a graph and an exons tab
rule all_index:
	input:
		os.path.join(config['whippet_bin'],config['run_name'],config['run_name'],".jls"),
		os.path.join(config['whippet_bin'],config['run_name'],config['run_name'],".exons.tab.gz")


rule merge_bam:
	input:
		bam_dir
	output:
		temp(os.path.join(config['data_output_path'],config['run_name'],"merged.bam"))
	shell:
		"""
		samtools merge bams
		"""
rule sort_bam:
	input:
		data_output_path + mergedbam
	output:
		temp(data_output_path + sorted_mergedbam)
	shell:
		"""
		samtools merge bamfiles
		"""
rule remove_dups:
	input:
		data_output_path + sorted_mergedbam
	output:
		temp(data_output_path + sorted_mergedbam_rm_dup)
	shell:
		"""
		samtools merge bamfiles
		"""
rule index_bam:
	input:
		data_output_path + sorted_mergedbam_rm_dup
	output:
		temp(data_output_path + sorted_mergedbam_rm_dup.bai)
	shell:
		"""
		samtools merge bamfiles
		"""
	
rule build_whippet_index:
	input:
		bam = data_output_path + sorted_mergedbam_rm_dup,
		bai = data_output_path + sorted_mergedbam_rm_dup.bai,
	output:
		whpt_graph = os.path.join(config['whippet_bin'],config['run_name'],config['run_name'],".jls"),
		whpt_exons = os.path.join(config['whippet_bin'],config['run_name'],config['run_name'],".exons.tab.gz")
	params:
		fasta = 
		gtf = 
		annotation_julia_indices = os.path.join(config['whippet_bin'],config['run_name'],config['run_name'])

	shell:
		"""
		config[julia] config[julia_bin]/whippet-index.jl --fasta {params.fasta} --bam {input.bam} --gtf {params.gtf} --index {params.annotation_julia_indices}
		"""


