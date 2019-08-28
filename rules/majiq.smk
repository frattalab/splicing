

import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"

rule buildMAJIQConfig:
	input:
		config['sample_csv_path']
	output:
		majiqConfig = config['majiq_top_level'] + config['run_name'] + "_majiqConfig.tsv"
	params:
		conditions_bams_parsed = parse_sample_csv_majiq(config['sample_csv_path'])
	# use config.yaml and samples.tsv to make MAJIQ file
	run:
		options = [
			"[info]",
			"readlen=" + str(config['read_len'])
			"samdir=" + config['read_len'],
			"genome=" + config['genome_refcode'],
			"strandness=" + config['strand_code'],
			"[experiments]"
		options += conditions_bams_parsed
		with open(test,"w") as outFile:
		    for opt in options:
		        outFile.write(opt + "\n")
rule majiq_build:
	input:
		majiq_config_file = config['majiq_config'] # then add the bams as in the input
	output:
		config['majiq_builder_output'] + "TDP_5.Aligned.sorted.out.splicegraph.sql."
	threads:
			16
	params:
		majiq_path = config['majiq_path'],
		gff3 = config['gff3'],
		majiq_builder_output = config['majiq_builder_output']
	shell:
		"""
		 {params.majiq_path} build {params.gff3} -c {input.majiq_config_file} -j {threads} -o {params.majiq_builder_output}
		"""
