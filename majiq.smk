

#todo impletent this
# # create MAJIQ config file
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
	output:
	threads:
			16
	shell:
		"""
		config['majiq_path'] build config['gff3'] -c config['majiq_config'] -j {threads} -o config['majiq_builder_output']
		"""
