import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"
include: "helpers.py"
# path to track and reference
TRACK   = config['gtf']
REF     = config['fasta']

#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))

print(SAMPLES)

CLASS1 = samples2.loc[samples.group == GROUPS[0]]['sample_name'].tolist()
CLASS2 = samples2.loc[samples.group == GROUPS[1]]['sample_name'].tolist()

# path to bam files
CLASS1_BAM = expand(config['bam_dir'] + '{sample}' + config['bam_suffix'] + '.bam', sample=CLASS1)
CLASS2_BAM = expand(config['bam_dir'] + '{sample}' + config['bam_suffix'] + '.bam', sample=CLASS2)



rule all:
    input:
        config['top_level_project_folder'] + 'cufflinks_merged/' + 'merged.gtf',
        config['top_level_project_folder'] + 'cufflinks_merged/' + 'comparison/all.stats',
        config['top_level_project_folder'] + 'diffexp/gene_exp.diff',
        config['top_level_project_folder'] + 'diffexp/isoform_exp.diff'


rule assembly:
    input:
        config['bam_dir'] + '{sample}' + '.bam'
    output:
        config['bam_dir'] + '{sample}/transcripts.gtf'
    wildcard_constraints:
        sample="|".join(SAMPLES)
    params:
        outdir = config['bam_dir'] + '{sample}',
        outmerged = config['bam_dir'] + 'cufflinks_merged/',
        strandness = get_cuff_strand(config['strand_code'])
    threads: 4
    shell:
        """
        mkdir -p {params.outmerged}
        mkdir -p {params.outdir}
        cufflinks --num-threads {threads} -o {params.outdir} --library-type {params.strandness} --frag-bias-correct {REF} {input}
        """


rule compose_merge:
    input:
        expand(config['bam_dir'] + '{sample}/transcripts.gtf', sample=SAMPLES)
    output:
        txt=config['bam_dir'] + 'cufflinks_merged/' + 'assemblies.txt'
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)


rule merge_assemblies:
    input:
        config['bam_dir'] + 'cufflinks_merged/' + 'assemblies.txt'
    output:
        config['bam_dir'] + 'cufflinks_merged/' + 'merged.gtf'
    params:
        dir=config['bam_dir'] + 'cufflinks_merged/'
    shell:
        """
        source /share/apps/source_files/python/python-2.7.16.source
        mkdir -p {params.dir}
        cuffmerge -o {params.dir} -s {REF} {input}
        """


rule compare_assemblies:
    input:
        config['bam_dir'] + 'cufflinks_merged/' + 'merged.gtf'
    output:
        config['bam_dir'] + 'cufflinks_merged/' + 'comparison/all.stats'
    params:
        dir=config['bam_dir'] + 'cufflinks_merged/' + 'comparison/'
    shell:
        """
        source /share/apps/source_files/python/python-2.7.16.source
        mkdir -p {params.dir}
        cuffcompare -o {params.dir}all -s {REF} -r {TRACK} {input}
        """

rule diffexp:
    input:
        class1=CLASS1_BAM,
        class2=CLASS2_BAM,
        gtf=config['bam_dir'] + 'cufflinks_merged/' + 'merged.gtf'
    output:
        config['bam_dir'] + 'diffexp/gene_exp.diff',
        config['bam_dir'] + 'diffexp/isoform_exp.diff'
    params:
        class1=",".join(CLASS1_BAM),
        class2=",".join(CLASS2_BAM),
        dir = config['bam_dir'] + 'diffexp/'
    threads: 8
    shell:
        """
        source /share/apps/source_files/python/python-2.7.16.source
        mkdir -p {params.dir}
        cuffdiff -o {params.dir} --num-threads {threads} {input.gtf} {params.class1} {params.class2}
        """
