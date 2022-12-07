# -*- coding: utf-8
"""
Created on 17:07 27/07/2018
Snakemake de novo transcriptomics and transcript abundance estimation
- Find the the fastq files used by DCC read alignments workflow
- compute de novo tx with StringTie [doi:10.1038/nprot.2016.095]
- use de novo annotation to compute transcript abundance with salmon [https://doi.org/10.1038/nmeth.4197]

If you use this workflow, please cite

Patro, R., et al. "Salmon provides fast and bias-aware quantification of transcript expression. Nat Meth. 2017; 14 (4): 417–9."
Pertea, Mihaela, et al. "Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown." Nature protocols 11.9 (2016): 1650.
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2019, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from itertools import groupby
import re


def extract_samples_replicates(samples, _pattern=re.compile("^(.+)_(.+)$")):
    """
    Extract pairs of condition and replicate name from sample files

    :param str _pattern: pattern to . Default uses {condition}_{replicate} template
    :param list samples:
    :return:
    :rtype: list
    """
    return list(zip(*[re.match(_pattern, x).groups() for x in samples]))


strand = {"fr-firststrand": "--fr", "fr-secondstrand": "--rf"}


workdir: config.get("path", ".")


container: "docker://tbrittoborges/stringtie:2.1.5"


cond, rep = extract_samples_replicates(config["samples"].keys())
name = config["samples"].keys()
raw_name = config["samples"].values()
sample_path = config["sample_path"]

d = {k: list(v) for k, v in groupby(sorted(zip(cond, rep)), key=lambda x: x[0])}

cond = set(cond)


include: "symlink.smk"


rule all:
    input:
        expand("mappings/{name}.bam", name=name),
        expand("stringtie/merged_bam/{group}.bam", group=cond),
        expand("stringtie/stringtie/{group}.gtf", group=cond),
        "stringtie/merged/merged.combined.gtf",


rule stringtie_merge_bam:
    input:
        lambda wc: ["mappings/{}_{}.bam".format(*x) for x in d[wc.group]],
    output:
        bam="stringtie/merged_bam/{group}.bam",
        bai="stringtie/merged_bam/{group}.bam.bai",
    threads: 10
    log:
        "logs/stringtie_merge_bam/{group}.log",
    wildcard_constraints:
        group="|".join(cond),
    envmodules:
        "samtools/1.12_deb10",
    shadow:
        "shallow"
    shell:
        "samtools merge {output.bam} {input} --threads {threads};samtools index {output.bam} {output.bai} 2> {log} "


rule stringtie_denovo_transcriptomics:
    input:
        "stringtie/merged_bam/{group}.bam",
    output:
        "stringtie/stringtie/{group}.gtf",
    params:
        strandness=strand.get(config.get("strandness", ""), ""),
        min_junct_coverage=config.get("min_junct_coverage", 3),
        min_isoform_proportion=config.get("min_isoform_proportion", 0.001),
        minimum_read_per_bp_coverage=config.get("minimum_read_per_bp_coverage", 3),
    wildcard_constraints:
        group="|".join(cond),
    log:
        "logs/stringtie_denovo_transcriptomics/{group}.log",
    envmodules:
        "stringtie",
    shadow:
        "shallow"
    shell:
        "stringtie {input} -o {output} -p {threads}  {params.strandness} -c {params.minimum_read_per_bp_coverage} -j {params.min_junct_coverage} -f {params.min_isoform_proportion} 2> {log} "


rule gffcompare:
    input:
        expand("stringtie/stringtie/{cond}.gtf", cond=cond),
    output:
        "stringtie/merged/merged.combined.gtf",
    log:
        "logs/gffcompare.log",
    params:
        gtf=config["ref"],
        out="stringtie/merged/merged",
    envmodules:
        "rnaseqtools",
    shadow:
        "shallow"
    shell:
        "gffcompare {input} -r {params.gtf} -R -V -o {params.out} 2>{log}"
