"""
Created on 17:07 27/02/2018
Snakemake workflow for JunctionSeq.

If you use JunctionSeq, please cite:

Hartley, Stephen W., and James C. Mullikin. "Detection and visualization of differential splicing in RNA-Seq data with
 JunctionSeq." Nucleic acids research 44.15 (2016): e127-e127.
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2020, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from itertools import groupby, chain


workdir: config.get("path", ".")


name = config["samples"].keys()
sample = config["samples"].values()
gtf_path = config["ref"]
comp_names = config["contrasts"].keys()
cond, rep = glob_wildcards("mappings/{cond}_{rep}.bam")
cond = set(cond)


container: "docker://tbrittoborges/junctionseq:1.16.0"


strandness = {"fr-secondstrand": "--stranded_fr_secondstrand", "fr-firststrand": "--stranded"}


rule all:
    input:
        "logs/",
        expand(
            "junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz",
            name=name,
        ),
        expand("junctionseq/{comparison}_decoder.tab", comparison=comp_names),
        "junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz",
        expand(
            "junctionseq/analysis/{comparison}_sigGenes.results.txt.gz",
            comparison=comp_names,
        ),


localrules:
    cat_decoder,
    create_decoder,


include: "symlink.smk"


rule junctionseq_qc:
    input:
        "mappings/{name}.bam",
    output:
        "junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz",
    params:
        gtf=config["ref"],
        output="junctionseq/rawCts/{name}/",
        strandness=strandness.get(config.get("strandness"), ""),
        max_read=config["read_len"],
        is_paired_end="--singleEnded" if config.get("is_single_end") == True else "",
    log:
        "logs/junctionseq_qc/{name}.log",
    envmodules:
        "java qorts ",
    shadow:
        "shallow"
    shell:
        "qorts QC {params.strandness} {params.is_paired_end} --maxReadLength {params.max_read} "
        "--runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon "
        "{input} {params.gtf} {params.output} 2> {log}"


rule junctionseq_create_decoder:
    log:
        "logs/junctionseq_create_decoder/{comparison}.log",
    output:
        "junctionseq/{comparison}_decoder.tab",
    run:
        ref, alt = wildcards[0].split("-vs-")
        with open(str(output), "w") as fou:
            fou.write("sample.ID\tgroup.ID\n")
            for n in name:
                if n.startswith(ref) or n.startswith(alt):
                    fou.write("{0}\t{1}\n".format(n, n.split("_")[0]))


rule junctionseq_cat_decoder:
    input:
        decoder=expand("junctionseq/{comparison}_decoder.tab", comparison=comp_names),
    output:
        "junctionseq/decoder.tab",
    log:
        "logs/junctionseq_cat_decoder.log",
    shadow:
        "shallow"
    shell:
        "awk 'FNR>1 || NR==1 ' {input.decoder} | awk '!x[$0]++' > {output} 2> {log}"


rule junctionseq_merge:
    input:
        decoder="junctionseq/decoder.tab",
        counts=expand(
            "junctionseq/rawCts/{name}/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz",
            name=name,
        ),
    output:
        "junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz",
    params:
        gtf=config["ref"],
        min_count=config.get("junctionseq_mincount", 6),
        strandness=strandness.get(config.get("strandness"), ""),
    envmodules:
        "java qorts",
    log:
        "logs/junctionseq_merge.log",
    shadow:
        "shallow"
    shell:
        "qorts mergeNovelSplices "
        "--minCount {params.min_count} "
        "{params.strandness} "
        "junctionseq/rawCts "
        "{input.decoder} "
        "{params.gtf} "
        "junctionseq/mergedOutput/ 2> {log}"


rule junctionseq_analysis:
    input:
        reads="junctionseq/mergedOutput/withNovel.forJunctionSeq.gff.gz",
        decoder="junctionseq/{comparison}_decoder.tab",
    output:
        "junctionseq/analysis/{comparison}_sigGenes.results.txt.gz",
    log:
        "logs/junctioseq_analysis/{comparison}.log",
    threads: 10
    resources:
        mem_mb=64000
    envmodules:
        "R/3.6.3_deb10 junctionseq/1.16.0_deb10",
    shadow:
        "shallow"
    shell:
        "junctionSeq.R {input.decoder} {output} {threads} 2> {log}"
