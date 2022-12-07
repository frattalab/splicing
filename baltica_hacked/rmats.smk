# -*- coding: utf-8
"""
Created on 17:07 27/07/2018
Snakemake workflow for rMATS

If you use rMATs, please cite

Shen, Shihao, et al. "rMATS: robust and flexible detection of differential 
alternative splicing from replicate RNA-Seq data." Proceedings of the 
National Academy of Sciences 111.51 (2014): E5593-E5601.
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2021, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from collections import defaultdict
import tempfile


container: "docker://tbrittoborges/rmats:latest"


workdir: config.get("path", ".")


contrasts = config["contrasts"]
keys = config["samples"].keys()
keys = [tuple(x.split("_")) for x in keys]
print("Hello beautiful")
temp_dir = tempfile.TemporaryDirectory()

d = defaultdict(list)
for x in keys:
    d[x[0]].append(x[1])

include: "symlink.smk"


localrules:
    create_rmats_input,


rule all:
    input:
        expand("rmats/{group}.txt", group=d.keys()),
        expand("rmats/{contrast}/", contrast=contrasts),


rule rmats_create_input:
    input:
        lambda wc: expand("mappings/{{group}}_{rep}.bam", rep=d.get(wc.group)),
    log:
        "logs/rmats/create_input_{group}.log",
    output:
        "rmats/{group}.txt",
    run:
        with open(str(output), "w") as fou:
            fou.write(",".join(input))


rule rmats_run:
    input:
        alt="rmats/{alt}.txt",
        ref="rmats/{ref}.txt",
    output:
        directory("rmats/{alt}-vs-{ref}/"),
    shadow:
        "shallow"
    log:
        "logs/rmats/run_{alt}-vs-{ref}.log",
    threads: 10
    envmodules:
        "rmats-turbo/4.1.1",
    params:
        gtf=config["ref"],
        is_paired="single" if config.get("is_single_end") else "paired",
        lib=config.get("strandness", "fr-unstranded"),
        read_len=config["read_len"],
        allow_clipping=config.get("rmats_allow_clipping", "--allow-clipping"),
        variable_read_length=config.get(
            "rmats_variable_read_length", "--variable-read-length"
        ),
        novel_ss=config.get("rmats_novel_ss", "--novelSS"),
        extra=config.get("rmats_extra", ""),
        tmp=os.path.join(temp_dir.name, "{alt}_vs_{ref}/"),
    shell:
        "rmats.py "
        "--b1 {input.alt} "
        "--b2 {input.ref} "
        "--gtf {params.gtf} "
        "--readLength {params.read_len} "
        "--nthread {threads} "
        "-t {params.is_paired} "
        "--libType {params.lib} "
        "{params.novel_ss} "
        "{params.allow_clipping} "
        "{params.variable_read_length} "
        "--od {output} "
        "--tmp {params.tmp} "
        "{params.extra}"
