# -*- coding: utf-8
"""
Created on 12:00 07/01/2020
Snakemake file for quality control of RNA-Sequencing dataset focused on junction quality.

Wang, Liguo, Shengqin Wang, and Wei Li. "RSeQC: quality control of RNA-seq experiments."
Bioinformatics 28.16 (2012): 2184-2185 and FastQC https://github.com/s-andrews/FastQC
"""

__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2020, Thiago Britto Borges"
__email__ = "tbrittoborges@uni-heidelberg.de"
__license__ = "MIT"

name = config["samples"].keys()
sample = config["samples"].values()


workdir: config.get("path", ".")


include: "symlink.smk"


container: "docker://tbrittoborges/qc:latest"


rule all:
    input:
        expand("qc/fastqc/{name}_fastqc.zip", name=name),
        "qc/ref.bed12",
        expand("qc/rseqc/{name}.inner_distance_plot.pdf", name=name),
        expand("qc/rseqc/{name}.GC_plot.pdf", name=name),
        expand("qc/rseqc/{name}.DupRate_plot.pdf", name=name),
        expand("qc/rseqc/{name}.splice_junction.pdf", name=name),
        expand("qc/rseqc/{name}.splice_events.pdf", name=name),
        expand("qc/rseqc/{name}.junctionSaturation_plot.pdf", name=name),
        expand("qc/rseqc/{name}.infer_experiment.txt", name=name),
        "qc/multiqc/multiqc_report.html",


# based on https://stackoverflow.com/a/50882104
rule fastqc:
    input:
        "mappings/{name}.bam",
    output:
        html="qc/fastqc/{name}_fastqc.html",
        zip="qc/fastqc/{name}_fastqc.zip",
    threads: 10
    params:
        prefix="qc/fastqc/",
    conda:
        "envs/qc.yml"
    envmodules:
        "fastqc",
    shadow:
        "shallow"
    log:
        "logs/fastqc_{name}.log",
    shell:
        "fastqc -t {threads} {input} -o {params.prefix} 2> {log}"


rule ref_annotation_gtf_to_bed:
    input:
        config["ref"],
    output:
        "qc/ref.bed12",
    conda:
        "envs/qc.yml"
    envmodules:
        "ucsc",
    shadow:
        "shallow"
    log:
        "logs/ref_annotation_gtf_to_bed.log",
    shell:
        # awk '{{if ($3 != "gene") print $0;}}' {input} | grep -v '^#' |
        """
        gtfToGenePred {input} /dev/stdout | genePredToBed /dev/stdin {output}
        """


rule rseqc_gene_body_coverage:
    input:
        bed="qc/ref.bed12",
    output:
        "qc/rseqc/geneBodyCoverage.curves.pdf",
    conda:
        "envs/qc.yml"
    envmodules:
        "rseqc",
    params:
        prefix="qc/rseqc/",
    shadow:
        "shallow"
    log:
        "logs/rseqc_gene_body_coverage.log",
    shell:
        "geneBody_coverage.py -r {input.bed} -i mappings/ -o {params.prefix} 2> {log}"


rule rseqc_inner_distance:
    input:
        bed="qc/ref.bed12",
        bam="mappings/{name}.bam",
    output:
        "qc/rseqc/{name}.inner_distance_plot.pdf",
    params:
        prefix="qc/rseqc/{name}",
    conda:
        "envs/qc.yml"
    envmodules:
        "rseqc",
    shadow:
        "shallow"
    log:
        "logs/rseqc_inner_distance_{name}.log",
    shell:
        "inner_distance.py -i {input.bam} -o {params.prefix} -r {input.bed} 2> {log}"


rule rseqc_read_gc:
    input:
        "mappings/{name}.bam",
    output:
        "qc/rseqc/{name}.GC_plot.pdf",
    params:
        prefix="qc/rseqc/{name}",
    envmodules:
        "rseqc",
    shadow:
        "shallow"
    log:
        "logs/rseqc_read_gc_{name}.log",
    shell:
        "read_GC.py -i {input} -o {params.prefix} 2> {log}"


### fix above output
rule rseqc_read_duplication:
    input:
        "mappings/{name}.bam",
    output:
        "qc/rseqc/{name}.DupRate_plot.pdf",
    params:
        prefix="qc/rseqc/{name}",
    conda:
        "envs/qc.yml"
    envmodules:
        "rseqc",
    shadow:
        "shallow"
    log:
        "logs/rseqc_read_duplication_{name}.log",
    shell:
        "read_duplication.py -i {input} -o {params.prefix} 2> {log}"


rule rseqc_junction_annotation:
    input:
        bed="qc/ref.bed12",
        bam="mappings/{name}.bam",
    output:
        "qc/rseqc/{name}.splice_junction.pdf",
        "qc/rseqc/{name}.splice_events.pdf",
    params:
        prefix="qc/rseqc/{name}",
    conda:
        "envs/qc.yml"
    envmodules:
        "rseqc",
    shadow:
        "shallow"
    log:
        "logs/rseqc_junction_annotation_{name}.log",
    shell:
        "junction_annotation.py -i {input.bam} -r {input.bed} -o {params.prefix} 2> {log}"


rule rseqc_junction_saturation:
    input:
        bed="qc/ref.bed12",
        bam="mappings/{name}.bam",
    output:
        "qc/rseqc/{name}.junctionSaturation_plot.pdf",
    params:
        prefix="qc/rseqc/{name}",
    envmodules:
        "rseqc",
    conda:
        "envs/qc.yml"
    shadow:
        "shallow"
    log:
        "logs/rseqc_junction_saturation_{name}.log",
    shell:
        "junction_saturation.py -i {input.bam} -r {input.bed} -o {params.prefix} 2> {log}"


rule rseqc_infer_experiment:
    input:
        bed="qc/ref.bed12",
        bam="mappings/{name}.bam",
    output:
        "qc/rseqc/{name}.infer_experiment.txt",
    conda:
        "envs/qc.yml"
    envmodules:
        "rseqc",
    shadow:
        "shallow"
    log:
        "logs/rseqc_infer_experiment_{name}.log",
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_bam_stat:
    input:
        "mappings/{name}.bam",
    output:
        "qc/rseqc/{name}.bam_stat.txt",
    envmodules:
        "rseqc",
    conda:
        "envs/qc.yml"
    shadow:
        "shallow"
    log:
        "logs/rseqc_bam_stat_{name}.log",
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_read_distribution:
    input:
        bam="mappings/{name}.bam",
        bed="qc/ref.bed12",
    envmodules:
        "rseqc",
    conda:
        "envs/qc.yml"
    output:
        "qc/rseqc/{name}.read_distribution.txt",
    shadow:
        "shallow"
    log:
        "logs/rseqc_read_distribution_{name}.log",
    shell:
        "read_distribution.py -i {input.bam} -r {input.bed} &> {output} 2> {log}"


rule multiqc:
    input:
        rules.all.input[:-1],
    envmodules:
        "multiqc",
    conda:
        "envs/qc.yml"
    output:
        "qc/multiqc/multiqc_report.html",
    shadow:
        "shallow"
    log:
        "logs/multiqc.log",
    shell:
        "multiqc -d --dirs-depth 1 qc/ -o qc/multiqc 2> {log} --force"
