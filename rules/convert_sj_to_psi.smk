import os
####GTF
gtf = "/SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v34.annotation.gtf"

####Folders and all the other stuff
####humans
out_spot = "normalized_annotated/"
# input_sj_folder = " /SAN/vyplab/alb_projects/data/sinai_splice_junctions/sinai_all_samples_renamed_sj_tabs/"
sj_suffix = ".SJ.out.tab"
####cell lines
input_sj_folder = "/SAN/vyplab/alb_projects/data/sinai_splice_junctions/all_bams_kds_linked/sj_files_only/"



# =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------

def get_single_psi_parsed_files_dasper(SAMPLES):
    """
    return a list of files that will exist
    """

    parsed_psi_files = [os.path.join(output_dir,x + "_normalized_annotated.csv") for x in SAMPLES]
    
    print(parsed_psi_files)
    return(parsed_psi_files)

output_dir = os.path.join(input_sj_folder,out_spot)
# print(bam_dir)
SAMPLES, = glob_wildcards(input_sj_folder + "{sample}" + sj_suffix)

rule all_normalize_annotate:
    input:
        expand(output_dir + "{sample}" + "_normalized_annotated.csv", sample = SAMPLES),
        os.path.join(output_dir, "combined_normalized_annotated.csv")


rule normalize_annotate:
    input:
        input_sj_folder + "{sample}" + sj_suffix
    output:
        output_dir + "{sample}" + "_normalized_annotated.csv"
    params:
        gtf = gtf,
        sample_name = "{sample}",
        output_folder = output_dir
    shell:
        """
        mkdir -p {output_dir}
        Rscript scripts/convert_sj_to_psi.R \
        --sample_name {params.sample_name} \
        --sample_file {input} \
        --gtf {params.gtf} \
        --output_folder {params.output_folder}
        """

rule squashed_normalize_annotate:
    input:
        all_parsed_csvs = get_single_psi_parsed_files_dasper(SAMPLES)
    output:
        os.path.join(output_dir, "combined_normalized_annotated.csv")
    params:
        dir_of_normed = output_dir
    shell:
        """
        Rscript scripts/combine_annotated_psi.R \
        --folder {params.dir_of_normed} \
        --out {output}
        """
