import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "../rules/helpers.py"

samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name']))
BASES, CONTRASTS = return_bases_and_contrasts()
ALLGROUP = list(set(BASES + CONTRASTS))
print(ALLGROUP)


SPECIES = config["species"]
GTF = config['gtf']

#make sure the output folder for STAR exists before running anything
bam_dir = get_output_dir(config["project_top_level"], config['bam_dir'])
bam_dir_temporary = get_output_dir(config["project_top_level"], 'merged_bams')

stringtie_outdir = get_output_dir(config["project_top_level"], config['stringtie_output'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])

both_output_dirs = [stringtie_outdir,scallop_outdir]


print(both_output_dirs)
rule allMerging:
    input:
        expand(os.path.join(bam_dir_temporary,"{grp}.bam"), grp = ALLGROUP),
        expand(os.path.join(bam_dir_temporary,"{grp}.bam.bai"),grp = ALLGROUP),
        expand(os.path.join(scallop_outdir,'{grp}' + ".gtf"),grp = ALLGROUP),
        expand(os.path.join(stringtie_outdir,'{grp}' + ".gtf"),grp = ALLGROUP),
        expand('{outputdir}{grp}' + ".gtf.map",outputdir =both_output_dirs, grp = ALLGROUP),

rule merge_bam_groups:
    input:
        group_bam_files = lambda wildcards: file_path_list(wildcards.grp,bam_dir,config['bam_suffix'] + '.bam')
    output:
        bam= os.path.join(bam_dir_temporary,"{grp}.bam"),
    wildcard_constraints:
        grp="|".join(ALLGROUP)
    threads: 10
    shell:
        """
        samtools merge {output.bam} {input.group_bam_files} --threads {threads}
        """

rule index_merged_group:
    input:
        bam = os.path.join(bam_dir_temporary,"{grp}.bam")
    output:
        bai= os.path.join(bam_dir_temporary,"{grp}.bam.bai")
    wildcard_constraints:
        grp="|".join(ALLGROUP)
    threads: 10
    shell:
        """
        samtools index {input.bam} {output.bai}
        """

rule scallop_per_group:
    input:
        bam= os.path.join(bam_dir_temporary,"{grp}.bam"),
        bai= os.path.join(bam_dir_temporary,"{grp}.bam.bai")
    output:
        os.path.join(scallop_outdir,'{grp}' + ".gtf")
    wildcard_constraints:
        grp="|".join(ALLGROUP)
    conda:
        "../envs/scallop.yml"
    params:
        scallop_path = 'scallop2',
        scallop_out_folder = scallop_outdir,
        scallop_extra_config = return_parsed_extra_params(config['scallop_extra_parameters']),
        libtype = config['scallop_strand']
    shell:
        """
        mkdir -p {params.scallop_out_folder}
        {params.scallop_path} \
        -i {input.bam} \
        -o {output} \
        --library_type {params.libtype} \
        {params.scallop_extra_config}
        """

rule stringtie_per_group:
    input:
        bam= os.path.join(bam_dir_temporary,"{grp}.bam"),
        bai= os.path.join(bam_dir_temporary,"{grp}.bam.bai")
    output:
        stringtie_outdir + "{grp}.gtf"
    wildcard_constraints:
        grp="|".join(ALLGROUP)
    conda:
        "../envs/stringtie.yml"
    shell:
        "stringtie {input.bam} -G {GTF} -o {output}"


rule compare_reference:
    input:
        "{outputdir}" + "{grp}.gtf"
    output:
        "{outputdir}{grp}.gtf.map"
    wildcard_constraints:
        grp="|".join(ALLGROUP)
    params:
        ref_gtf = GTF,
        gffcompare = config['gffcompare'],
        prefix = os.path.join('{outputdir}', '{grp}')
    shell:
        """
        {params.gffcompare} -o gffall -r {params.ref_gtf} {input}
        """

# # rule fetch_unique:
# #     input:
# #         sample_tmap = os.path.join(scallop_outdir,"scallop_merged.gtf"),
# #         sample_gtf = os.path.join(scallop_outdir, "gffall.scallop_merged.gtf.tmap")
# #     output:
# #         os.path.join(scallop_outdir, "scallop_merged.unique.gtf")
# #     params:
# #         ref_gtf = GTF,
# #         gtfcuff = config['gtfcuff']
# #     shell:
# #         """
# #         {params.gtfcuff} puniq {input.sample_tmap} {input.sample_gtf} {params.ref_gtf} {output}
# #         """
