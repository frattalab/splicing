
rule merge_bam:
	input:
		all the bam file
	output:
		temp(data_output_path + mergedbam)
	shell:
		"""
		smatools merge bamfiles
		"""
rule sort_bam:
	input:
		data_output_path + mergedbam
	output:
		temp(data_output_path + sorted_mergedbam)
	shell:
		"""
		smatools merge bamfiles
		"""
rule remove_dups:
	input:
		data_output_path + sorted_mergedbam
	output:
		temp(data_output_path + sorted_mergedbam_rm_dup)
	shell:
		"""
		smatools merge bamfiles
		"""
rule index_bam:
	input:
		data_output_path + sorted_mergedbam_rm_dup
	output:
		temp(data_output_path + sorted_mergedbam_rm_dup.bai)
	shell:
		"""
		smatools merge bamfiles
		"""
	
rule build_whippet_index:
	input:
		bam = data_output_path + sorted_mergedbam_rm_dup,
		bai = data_output_path + sorted_mergedbam_rm_dup.bai,
	output:
		whippet bin + run_name_index,
		temp
	params:
		fasta = 
		gtf = 

	shell:


