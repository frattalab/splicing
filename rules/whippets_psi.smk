import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"


#include: "rules/fastp.smk"
#RULE ORDER DIRECTIVE
#if paired end, use the paired end rule to run, if single end use the single end rule to run
# if config['end_type'] == "pe":
# 	ruleorder: run_star_pe > run_star_se
# else:
# 	ruleorder: run_star_se > run_star_pe


SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

###make the whippet output final FOLDERR
output_dir = os.path.join(config['top_level_project_folder'],config['whippet_output_path'])
#make sure the output folder for Whippets exists before running anything
os.system("mkdir -p {0}".format(output_dir))

rule all_whip:
    input:
        expand(output_dir + "{sample}.psi.gz", sample = SAMPLES)
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)

if config['endtype'] == "se":
    rule whippet_psi:
        input:
            fq = config['fastq_files'] + "{sample}_1.merged.fastq.gz",
            index = os.path.join(whippet_index_path + "{sample}" + ".jls")
    rule whippet_psi:
        input:
            fast1 = fastq_dir  + "{sample}_1.merged.fastq.gz",
            fast2 = fastq_dir  + "{sample}_2.merged.fastq.gz",
        output:
            output_dir + "{sample}/" + "quant.sf"
        params:
            output_dir + "{sample}"
        threads: 2
        shell:
            """
            /share/apps/julia-0.6.3/bin/julia /SAN/vyplab/alb_projects/tools/julia_pkgdir/v0.6/Whippet/bin/whippet-quant.jl \
            {input.fast1} \
            {input.fast2} \
            -x /home/annbrown/data/ward_bams/whippets/scallop_merged_first.jls\
            --out {params}
             """

#$ -l tmem=16G
#$ -l h_vmem=16G
#$ -l h_rt=72:0:0

#$ -j y
#$ -N xwhip
