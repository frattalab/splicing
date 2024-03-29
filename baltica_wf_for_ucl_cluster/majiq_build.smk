# -*- coding: utf-8
"""
Created on 17:07 29/02/2018 2018
Snakemake workflow for Majiq.

If you use Majiq, please cite

Vaquero-Garcia, Jorge, et al. "A new view of transcriptome complexity and
regulation through the lens of local splicing variations." elife 5 (2016):
e11752.
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2020, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"


def comparison(wc, index, mapping):
    condition = wc.contrast.split("-vs-")[index]
    return expand("majiq/{name}.majiq", name=mapping[condition])


def natural_sort_key(s, _nsre=re.compile("([0-9]+)")):
    return [int(text) if text.isdigit() else text.lower() for text in _nsre.split(s)]

strand = {"fr-firststrand": "reverse", "fr-secondstrand": "forward"}
strand = strand.get(config.get("strandness"), "None")
def create_ini(output):
    lines = [
        "[info]",
        "readlen={}".format(config["read_len"]),
        "bamdirs={}".format("mappings/"),
        "genome={}".format(config.get("assembly", "")),
        "genome_path={}".format(config["ref_fa"]),
        "strandness={}".format(strand),
        "[experiments]",
    ]

    lines.extend(["{}={}".format(k, ",".join(v)) for k, v in mapping.items()])

    with open(str(output), "w") as ini:
        ini.writelines("\n".join(lines))


container: "docker://tbrittoborges/majiq:2.2"


workdir: config.get("path", ".")



name = config["samples"].keys()
sample = config["samples"].values()
sample_path = config["sample_path"]
contrasts = config["contrasts"]
conditions = sorted(set([x.split("_")[0] for x in name]), key=natural_sort_key)

mapping = {c: [x for x in name if x[: x.index("_")] == c] for c in conditions}


localrules:
    symlink,
    majiq_create_ini,
    majiq_gtf_to_gff


include: "symlink.smk"


rule all:
    input:
        "logs/",
        expand("mappings/{name}.bam", name=name),
        "majiq/build.ini",
        expand("majiq/{name}.majiq", name=name)


rule majiq_create_ini:
    input:
        expand("mappings/{name}.bam", name=name),
    output:
        "majiq/build.ini",
    log:
        "logs/majiq_create_ini.log",
    run:
        create_ini(output)


rule majiq_gtf_to_gff:
    input: config["gff3"]
    output: "majiq/ref.gff"
    shell:
        "ln -s {input} {output}"


rule majiq_build:
    input:
        ini="majiq/build.ini",
        ref="majiq/ref.gff",
    output:
        expand("majiq/{name}.majiq", name=name),
        "majiq/splicegraph.sql",
    threads: len(conditions)
    log:
        "logs/majiq_build.log",
    shell:
        " majiq build --conf {input.ini} --nproc {threads} --output majiq/ {input.ref}"

