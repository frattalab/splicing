

import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"

# rule buildMAJIQConfig:
# 	output:
# 		majiqConfig = outFolder + dataCode + "_majiqConfig.tsv"
# 	# use config.yaml and samples.tsv to make MAJIQ file
# 	run:
# 	# import yaml
# 		# fix bamSuffix - default is ".bam" but often have other strings
# 		# ideally keep same as for leafcutter pipeline - less to worry about
# 		options = [
# 			"[info]",
# 			"readlen=" + str(config['readLen'])
# 			"samdir=" + config['readLen'],
# 			"genome=" + config['refcode'],
# 			"strandness=" + config['strandCode'],
# 			"[experiments]",
# 			refCondition + "=" + (bamSuffix + ",").join(refSamples) + bamSuffix,
# 			altCondition + "=" + (bamSuffix + ",").join(altSamples) + bamSuffix
# 			]
# 		with open(output.majiqConfig,"w") as outFile:
# 			for opt in options:
# 				outFile.write(opt + "\n")
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
