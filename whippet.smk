

import os,re
#folder where the fastq files are
fastq_files = "/Users/annaleigh/Documents/data/4su_bams/fastq/"
#lazy way of pulling out all the sample names in the fastq fiels
SAMPLES = [re.sub(".merged.fastq.gz", "",smp) for smp in os.listdir(fastq_files)]
#eventual output folder for the whippet analysis
data_output_path = "/SAN/vyplab/alb_projects/data/muscle/analysis/whippets/"
#full path to your julia install
julia_path = "/Applications/Julia-0.6.app/Contents/Resources/julia/bin/julia"
#the path for the whippset-quant jl script
whippet_bin = "/home/annbrown/alb_projects/tools/julia_pkgdir/v0.6/Whippet/bin/"
#the path where your samtools is installed
samtools_path = "/share/apps/genomics/samtools-1.9/bin/samtools"
#the path to the julia index (we're going to make a rule later to define this)
INDEX = "/Users/annaleigh/.julia/v0.6/Whippet/index/with_merged_bams.jls"
#first rule of a snakemake file always defines the final desired output. Eventually this file should produced the psi file
#and the sorted bam 
print(SAMPLES)

#our final output
rule all_whippet:
	input:
		expand(data_output_path + "{sample}.sorted.bam.bai", sample = SAMPLES)
#whippet's need to have an index for a project in order to run properly
#in order to include novel junctions, we can build it using the aligned bam files
#



#first we'll run the psi runn
rule whippet_psi:
	input:
		"/Users/annaleigh/Documents/data/4su_bams/fastq/{sample}.merged.fastq.gz"
		 #input is the merged fastq files for all the samples
	output:
		psi_file = data_output_path + "{sample}.psi.gz",
		jnc_file = data_output_path + "{sample}.jnc.gz",
		#i believe this marks the output file as temporary so that it is deleted one every consuming job has been executed

		sam_file = temp(data_output_path + "{sample}.sam")
	params:
		julia = julia_path,
		whippet_quant = whippet_quant,
		index = INDEX,
		output_path = lambda wildcards: data_output_path + wildcards.sample # this is a really simple rule just run the whippet quant with the selected index the name of the output path and put it as a sam file
	shell:
		"""
		{params.julia} {params.whippet_quant} {input} --index {params.index} --o {params.output_path} --sam > {output.sam_file}
		"""
#now a quick set of rule to turn all those sam fiels to bam files with samtools; very standard stuff here
rule sam_to_bam:
	input:
		sam_file = data_output_path + "{sample}.sam"
	output:
	#i believe this marks the output file as temporary so that it is deleted one every consuming job has been executed
		temp(data_output_path + "{sample}.bam")
	params:
		samtools = samtools_path
	shell:
		"""
		{params.samtools} view -b {input.sam_file} > {output} 
		"""
rule sort_bams:
	input:
		data_output_path + "{sample}.bam"
	output:
		data_output_path + "{sample}.sorted.bam"
	params:
		samtools = samtools_path,
		old_sam = data_output_path + "{sample}.sam"
	shell:
		"""
		rm {params.old_sam}
		{params.samtools} sort {input} > {output}
		"""

rule index_sort_bams:
	input:
		data_output_path + "{sample}.sorted.bam"
	output:
		data_output_path + "{sample}.sorted.bam.bai"
	params:
		samtools = samtools_path,
		unsorted_bam = data_output_path + "{sample}.bam"
	shell:
		"""
		rm {params.unsorted_bam}
		{params.samtools} index {input} > {output}
		"""










