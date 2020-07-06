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


SAMPLES = pd.read_csv(config["sample_csv_path"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)
SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

###make the whippet output final FOLDERR
output_dir = os.path.join(config['top_level_project_folder'],config['whippet_output_path'])
whippet_psi_path = os.path.join(config['top_level_project_folder'], config['whippet_output_path'], "psi/")

#make sure the output folder for Whippets exists before running anything
os.system("mkdir -p {0}".format(output_dir))
os.system("mkdir -p {0}".format(whippet_psi_path))

rule all_whip:
    input:
        expand(whippet_psi_path + "{sample}.psi.gz", sample = SAMPLE_NAMES),
        expand(whippet_psi_path + "{sample}.jnc.gz", sample = SAMPLE_NAMES)

    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)

if config['endtype'] == "pe":
    rule whippet_psi:
        input:
            #fq = config['fastq_files'] + "{sample}_1.merged.fastq.gz",
            index = os.path.join(whippet_index_path, "{sample}" + ".jls"),
            fast1 = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = True),
            fast2 = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = False),

        output:
            psi = os.path.join(whippet_psi_path,"{sample}.psi.gz"),
            jnc_file = os.path.join(whippet_psi_path,"{sample}.jnc.gz"),
            sam_file = temp(os.path.join(whippet_psi_path,"{sample}.sam"))

        params:
            julia = config['julia'],
            quant_script = os.path.join(config['whippet_bin'],"whippet-quant.jl"),
            output_prefix = os.path.join(whippet_psi_path,"{sample}")

        threads: 2

        shell:
            """
            export JULIA_PKGDIR=/SAN/vyplab/alb_projects/tools/julia_pkgdir/v0.6/
            {params.julia} {params.quant_script} \
            {input.fast1} \
            {input.fast2} \
            -x {input.index} \
            -o {params.output_prefix}
            """

elif config['endtype'] == "se":
    rule whippet_psi:
        input:
            fastq = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = True),
            index = os.path.join(whippet_index_path, "{sample}" + ".jls"),

        output:
            psi = os.path.join(whippet_psi_path,"{sample}.psi.gz"),
            jnc_file = os.path.join(whippet_psi_path,"{sample}.jnc.gz"),
            sam_file = temp(os.path.join(whippet_psi_path,"{sample}.sam"))

        params:
            julia = config['julia'],
            quant_script = os.path.join(config['whippet_bin'],"whippet-quant.jl"),
            output_prefix = os.path.join(whippet_psi_path,"{sample}")

        threads: 2

        shell:
            """
            export JULIA_PKGDIR=/SAN/vyplab/alb_projects/tools/julia_pkgdir/v0.6/
            {params.julia} {params.quant_script} \
            {input.fastq} \
            -x {input.index} \
            --out {params.output_prefix}
            """
#$ -l tmem=16G
#$ -l h_vmem=16G
#$ -l h_rt=72:0:0

#$ -j y
#$ -N xwhip
