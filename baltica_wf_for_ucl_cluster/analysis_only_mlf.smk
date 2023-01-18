# -*- coding: utf-8
"""
Created on 20/05/2020
Snakemake workflow for analysis within Baltica.

usage:
    snakemake -s majiq.smk --configfile config.yaml -j 10
"""


workdir: config.get("path", ".")


name = config["samples"].keys()
contrasts = config["contrasts"]
project_title = config.get("project_title", "").replace(' ', '_')

container: "docker://tbrittoborges/baltica_analysis:1.0"


import sys
import pathlib

fastqc_file = "qc/multiqc/multiqc_data/multiqc_fastqc.txt"


rule all:
    input:
        "results/SJ_annotated.csv",
        "results/SJ_annotated_assigned.csv",
        expand(
            "results/baltica_report{project_title}.html", 
            project_title="_" + project_title),
        expand(
            "results/baltica_table{project_title}.xlsx", 
            project_title="_" + project_title),




rule parse_majiq:
    input:
        expand("majiq/voila/{contrast}_voila.tsv", contrast=config["contrasts"].keys()),
    output:
        "majiq/majiq_junctions.csv",
    envmodules:
        "R/4.0.5_deb10",
    log:
        "logs/parse_majiq.log",
    params:
        cutoff=1.1,
    script:
        "parse_majiq_output.R"


rule parse_leafcutter:
    input:
        expand(
            "leafcutter/{contrast}/{contrast}_cluster_significance.txt",
            contrast=contrasts.keys(),
        ),
    output:
        "leafcutter/leafcutter_junctions.csv",
    envmodules:
        "R/4.0.5_deb10",
    log:
        "logs/parse_leafcutter.log",
    params:
        cutoff=1.1,
    script:
        "parse_leafcutter_output.R"


rule annotate:
    input:
        expand(
            "{method}/{method}_junctions.csv",
            method=["majiq", "leafcutter"],
        ),
        ref="stringtie/merged/merged.combined.gtf",
    params:
        ref=config.get("ref"),
        orthogonal_result=config.get('orthogonal_result'),
        unstranded=False if 'strandness' in config else True,
    envmodules:
        "R/4.0.5_deb10",
    log:
        "logs/baltica_annotate.log",
    output:
        "results/SJ_annotated.csv",
    script:
        "annotate_SJ.R"
