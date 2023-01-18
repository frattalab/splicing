# -*- coding: utf-8
"""
Created on 17:07 29/02/2018 2018
Snakemake workflow for LeafCutter.

If you use LeafCutter, please cite:

Li, Yang I., et al. "Annotation-free quantification of RNA splicing using LeafCutter." 
Nature genetics 50.1 (2018): 151-158.

..Usage:
    snakemake -s leafcutter.smk --configfile config.yml
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2020, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from pathlib import Path
from subprocess import CalledProcessError

def basename(path, suffix=None):
    if suffix:
        return str(Path(path).with_suffix(suffix).name)
    return str(Path(path).name)


container: "docker://tbrittoborges/leafcutter:latest"


workdir: config.get("path", ".")


name = config["samples"].keys()
sample = config["samples"].values()
gtf_path = config["ref"]
conditions = [x.split("_")[0] for x in name]
sample_path = config["sample_path"]
comp_names = config["contrasts"].keys()

strand = {"fr-secondstrand": 2, "fr-firststrand": 1}


localrules:
    all,
    leafcutter_concatenate,
    symlink,


include: "symlink.smk"


rule all:
    input:
        "logs/",
        expand("mappings/{name}.bam", name=name),
        expand(
            "leafcutter/{comp_names}/{comp_names}_cluster_significance.txt",
            comp_names=comp_names,
        ),
        expand("leafcutter/{name}.junc", name=name),


rule leafcutter_bam2junc:
    # Details https://regtools.readthedocs.io/en/latest/commands/junctions-extract/
    input:
        "mappings/{name}.bam",
    output:
        "leafcutter/{name}.junc",
    envmodules:
        "regtools/0.6",
    log:
        "logs/leafcutter/leafcutter_bam2junc/{name}.log",
    params:
        # defaults as suggested by leafcutter
        minimum_anchor_length=config.get("leafcutter_minimum_anchor_length", 8),
        minimum_intron_size=config.get("leafcutter_minimum_intron_size", 50),
        maximum_intron_size=config.get("leafcutter_maximum_intron_size", 5e5),
        # Strand specificity of RNA library preparation, 
        # where 0 = unstranded/XS, 
        # 1 = first-strand/RF, 
        # 2 = second-strand/FR
        # If your alignments contain XS tags,
        # these will be used in the "unstranded" mode.
        strand_specificity=strand.get(config.get("strandness", 2), 0),
    envmodules:
        "regtools",
    shadow:
        "shallow"
    shell:
        """
        regtools junctions extract \
        -a {params.minimum_anchor_length} \
        -i {params.minimum_intron_size} \
        -I {params.maximum_intron_size} \
        -s {params.strand_specificity} \
        {input} > {output} 2> {log}
        """


rule leafcutter_concatenate:
    input:
        expand(rules.leafcutter_bam2junc.output, name=name),
    log:
        "logs/leafcutter/concatenate_{comp_names}.log",
    output:
        junc="leafcutter/{comp_names}/juncfiles.txt",
        test="leafcutter/{comp_names}/diff_introns.txt",
    run:
        comp = output.junc.split("/")[1]
        cond_a, cond_b = comp.split("-vs-")
        with open(output.junc, "w") as file_out:
            work_path = Path(config.get("path", "."))
            for n in name:
                if n.startswith(cond_a) or n.startswith(cond_b):
                    file_out.write(str(work_path / "leafcutter/{}.junc\n".format(n)))

        with open(output.test, "w") as file_out:
            for n in name:
                if n.startswith(cond_a):
                    file_out.write("{} {}\n".format(n, cond_a))

                elif n.startswith(cond_b):
                    file_out.write("{} {}\n".format(n, cond_b))


# step 2
# m= min numb reads per cluster, l = intron length
rule leafcutter_intron_clustering:
    input:
        rules.leafcutter_concatenate.output.junc,
    params:
        m=config.get("leafcutter_min_cluster_reads", 50),
        l=config.get("leafcutter_max_intron_length", 500000),
        prefix="leafcutter/{comp_names}/{comp_names}",
        n="{comp_names}",
    output:
        "leafcutter/{comp_names}/{comp_names}_perind_numers.counts.gz",
    log:
        "logs/leafcutter/leafcutter_intron_clustering/{comp_names}.log",
    envmodules:
        "python3/3.6.13_deb10",
    shadow:
        "shallow"
    shell:
        """
        leafcutter_cluster_regtools_py3.py \
        -j {input} -m {params.m} \
        -o {params.prefix} \
        -l {params.l}  2> {log}
        """


# step 3.1
rule leafcutter_gtf_to_exon:
    input:
        gtf_path,
    output:
        a="leafcutter/{}".format(basename(gtf_path, suffix=".gz")),
        b="leafcutter/exons.gtf.gz",
    envmodules:
        "R/4.0.5_deb10 leafcutter/0.2.7_deb10",
    log:
        "logs/leafcutter/leafcutter_gtf_to_exon.log",
    shadow:
        "shallow"
    shell:
        """
        gzip -c {input} > {output.a}
        gtf_to_exons.R {output.a} {output.b} 2> {log}
        """


rule leafcutter_differential_splicing:
    input:
        a=rules.leafcutter_gtf_to_exon.output.b,
        b="leafcutter/{comp_names}/{comp_names}_perind_numers.counts.gz",
        c="leafcutter/{comp_names}/diff_introns.txt",
    output:
        "leafcutter/{comp_names}/{comp_names}_cluster_significance.txt",
    params:
        min_samples_per_group=config.get("leafcutter_min_samples_per_group", 3),
        min_samples_per_intron=config.get("leafcutter_min_samples_per_intron", 5),
        min_coverage=config.get("leafcutter_min_coverage", 20),
        prefix="leafcutter/{comp_names}/{comp_names}",
    threads: 10
    conda:
        "../envs/leafcutter.yml"
    envmodules:
        "R/4.0.5_deb10 leafcutter/0.2.7_deb10",
    log:
        "logs/leafcutter/leafcutter_differential_splicing/{comp_names}.log",
    shadow:
        "shallow"
    shell:
        """
        leafcutter_ds_pair.R --exon_file={input.a} \
        --min_coverage {params.min_coverage}  \
        {input.b} {input.c} --num_threads {threads} --output_prefix={params.prefix} \
        -i {params.min_samples_per_intron} -g {params.min_samples_per_group} 2> {log}
        """

onsuccess:
    try:
        shell("rm *.sorted.gz")
    except CalledProcessError:
        pass