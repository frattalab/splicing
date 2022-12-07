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


rule parse_junctionseq:
    input:
        expand(
            "junctionseq/analysis/{contrast}_sigGenes.results.txt.gz",
            contrast=contrasts.keys(),
        ),
    output:
        "junctionseq/junctionseq_junctions.csv",
    envmodules:
        "R/4.0.5_deb10",
    log:
        "logs/parse_junctionseq.log",
    params:
        cutoff=1.1,
    script:
        "parse_junctionseq_output.R"


rule parse_rmats:
    input:
        expand(
            "rmats/{contrast}/{st}.MATS.JC.txt",
            contrast=contrasts.keys(),
            st=["A3SS", "A5SS", "RI", "MXE", "SE"],
        ),
    output:
        "rmats/rmats_junctions.csv",
    envmodules:
        "R/4.0.5_deb10",
    log:
        "logs/parse_rmats.log",
    params:
        cutoff=1.1,
    log:
        "logs/parse_rmats.log",
    script:
        "parse_rmats_output.R"


rule annotate:
    input:
        expand(
            "{method}/{method}_junctions.csv",
            method=["majiq", "leafcutter", "junctionseq", "rmats"],
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


rule assign_AS_type:
    input:
        "results/SJ_annotated.csv",
        ref="stringtie/merged/merged.combined.gtf",
    envmodules:
        "R/4.0.5_deb10",
    output:
        "results/SJ_annotated_assigned.csv",
    log:
        "logs/assign_AS_type.log",
    script:
        "assign_AS_type.R"


rule simplify:
    input:
        "results/SJ_annotated_assigned.csv",
    output:
        expand(
            "results/baltica_table{project_title}.xlsx", 
            project_title="_" + project_title),
    envmodules:
        "R/4.0.5_deb10",
    log:
        "logs/simplify.log",
    script:
        "simplify.R"


rule baltica_report:
    input:
        "results/SJ_annotated.csv",
        "leafcutter/leafcutter_junctions.csv"
    params:
        fastqc=fastqc_file if os.path.isfile(fastqc_file) else None,
        config=config["config_path"],
        # TODO streamline start_sj_file
        # star_sj=config.get("star_sj_file"),
        doc_title=config.get("project_title", ""),
        doc_authors=config.get("project_authors", ""),
    output:
        expand(
            "results/baltica_report{project_title}.html", 
            project_title="_" + project_title),
    log:
        "logs/baltica_report.log",    
    script:
        "baltica_report.Rmd"
