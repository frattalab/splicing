
import pandas as pd
import os
import subprocess
import yaml

configfile: "config/config.yaml"
include: "helpers.py"

if config['species'] == "mouse":
    sg_seq_output_folder = "/SAN/vyplab/vyplab_reference_genomes/SGSeq/mouse/"
elif config['species'] == "human":
    sg_seq_output_folder = "/SAN/vyplab/vyplab_reference_genomes/SGSeq/human/"
else:
    print("species currently set up for either 'mouse' or 'human'")

rule NAME:
    input:
        config['gtf']
    params:
        sg_seq_output
    threads:
        4
    output:
        sg_seq_output_folder = ""
    script:
        "scripts/sg_seq_create_anno.R"
