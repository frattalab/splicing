

import pandas as pd
import os
import subprocess
import yaml

configfile: "config/config.yaml"
include: "helpers.py"
samples = pd.read_csv(config['sample_csv_path'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name']))

whippet_index_path = config['whippet_output_path'] + "index/"
whippet_psi_path = config['whippet_output_path'] + "psi/"
# What a horrible way to make these directories, AL be ashamed XD
cmd = "mkdir -p {0}".format(whippet_index_path)
os.system(cmd)
cmd = "mkdir -p {0}".format(whippet_psi_path)
os.system(cmd)


# had conversation with Ben Whathisface at RNA uk, suggested trying to
# build each bam their own index, and using the PSI for that
# that does mean that you can't use the delte PSI, but let's give her a go
os.system("export JULIA_PKGDIR=/SAN/vyplab/alb_projects/tools/julia_pkgdir/v0.6/")

print(SAMPLE_NAMES)
print(len(SAMPLE_NAMES))
#our final output
rule all_whip:
    input:
        expand(whippet_index_path + "{sample}" + ".jls",sample = SAMPLE_NAMES),
        expand(whippet_psi_path + "{sample}.sorted.bam.bai", sample = SAMPLE_NAMES)
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
        # expand(os.path.join(config['whippet_bin'],"index", config['run_name'] + "gencode" + ".exons.tab.gz"))

rule build_whippet_index:
    input:
        bam = os.path.join(config['bam_dir'], "{sample}" + config['bam_suffix'] + ".bam"),
        bai = os.path.join(config['bam_dir'],  "{sample}" + config['bam_suffix'] + ".bam.bai")
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    output:
        whpt_graph = os.path.join(whippet_index_path + "{sample}" + ".jls")
    params:
        fasta = config['fasta'],
        gtf = config['gtf'],
        extra_whippet_params = return_parsed_extra_params(config['whippet_extra_parameters']),
        whippet_ind_path = os.path.join(config['whippet_bin'],"whippet-index.jl"),
        annotation_julia_indices = os.path.join(config['whippet_bin'],"index",config['run_name'])
    shell:
        """
        export JULIA_PKGDIR=/SAN/vyplab/alb_projects/tools/julia_pkgdir/v0.6/
        {config[julia]} {params.whippet_ind_path} --fasta {params.fasta} {params.extra_whippet_params} --bam {input.bam} --gtf {params.gtf} --index {output.whpt_graph}
        """
#
#first we'll run the psi run
if config['endtype'] == "se":
    rule whippet_psi:
        input:
            config['fastq_files'] + "{sample}_1.merged.fastq.gz"
        output:
            psi_file = whippet_psi_path + "{sample}.psi.gz",
            jnc_file = whippet_psi_path + "{sample}.jnc.gz",
            sam_file = temp(whippet_psi_path + "{sample}.sam")
        params:
            julia = config['julia'],
            whippet_quant = config['whippet_bin'] + "whippet-quant.jl",
            index = os.path.join(whippet_index_path + "{sample}" + ".jls"),
            output_path = whippet_psi_path
        shell:
            """
            export JULIA_PKGDIR=/SAN/vyplab/alb_projects/tools/julia_pkgdir/v0.6/
            {params.julia} {params.whippet_quant} {input} --index {params.index} --o {params.output_path} --sam > {output.sam_file}
            """
elif config['endtype'] == "pe":
        rule whippet_psi:
            input:
                fwd = config['fastq_files'] + "{sample}_1.merged.fastq.gz",
                rev = config['fastq_files'] + "{sample}_2.merged.fastq.gz"
            output:
                psi_file = whippet_psi_path + "{sample}.psi.gz",
                jnc_file = whippet_psi_path + "{sample}.jnc.gz",
                sam_file = temp(whippet_psi_path + "{sample}.sam")
            params:
                julia = config['julia'],
                whippet_quant = config['whippet_bin'] + "whippet-quant.jl",
                index = os.path.join(whippet_index_path + "{sample}" + ".jls"),
                output_path = whippet_psi_path
            shell:
                """
                export JULIA_PKGDIR=/SAN/vyplab/alb_projects/tools/julia_pkgdir/v0.6/
                {params.julia} {params.whippet_quant} {input} --index {params.index} --o {params.output_path} --sam > {output.sam_file}
                """
else:
    print("End type either 'se' or 'pe'")
#now a quick set of rule to turn all those sam fiels to bam files with samtools; very standard stuff here
rule sam_to_bam:
    input:
        sam_file = whippet_psi_path + "{sample}.sam"
    output:
    #i believe this marks the output file as temporary so that it is deleted one every consuming job has been executed
        temp(whippet_psi_path + "{sample}.bam")
    params:
        samtools = config['samtools_path']
    shell:
        """
        {params.samtools} view -b {input.sam_file} > {output}
        """
rule sort_bams:
    input:
        whippet_psi_path + "{sample}.bam"
    output:
        whippet_psi_path + "{sample}.sorted.bam"
    params:
        samtools = config['samtools_path'],
        old_sam = whippet_psi_path + "{sample}.sam"
    shell:
        """
        rm {params.old_sam}
        {params.samtools} sort {input} > {output}
        """

rule index_sort_bams:
    input:
        whippet_psi_path + "{sample}.sorted.bam"
    output:
        whippet_psi_path + "{sample}.sorted.bam.bai"
    params:
        samtools = config['samtools_path'],
        unsorted_bam = whippet_psi_path + "{sample}.bam"
    shell:
        """
        rm {params.unsorted_bam}
        {params.samtools} index {input} > {output}
        """
