import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"

#final output of the index run is a graph and an exons tab
rule all_index:
    input:
        os.path.join(config['whippet_indicies'],config['whippet_index_name'] + ".jls"),
        os.path.join(config['whippet_indicies'],config['whippet_index_name'] + ".exons.tab.gz")

rule build_whippet_index:
    input:
        fasta = config['fasta'],
        gtf = config['gtf'],
    output:
        whpt_graph = os.path.join(config['whippet_indicies'],config['whippet_index_name'] + ".jls"),
        whpt_exons = os.path.join(config['whippet_indicies'],config['whippet_index_name'] + ".exons.tab.gz")
    params:
        whippet_ind_path = os.path.join(config['whippet_bin'],"whippet-index.jl"),
        annotation_julia_indices = os.path.join(config['whippet_indicies'],config['whippet_index_name'])
    shell:
        """
        {config[julia]} {params.whippet_ind_path} \
        --fasta {params.fasta} \
        --gtf {params.gtf} \
        --index {params.annotation_julia_indices}
        """
