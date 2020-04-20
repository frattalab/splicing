import os
import pandas as pd
import yaml

def load_comparisons():
    comparisons = "config/comparisons.yaml"
    with open(comparisons, 'r') as stream:
        try:
            compare_dict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return(compare_dict)

def check_key(dict, key):
     """
     simple funciton to check if a key is in a dictionary, if it's not returns false and if it is returns the value
     """
     if key in dict:
         return(dict[key])
     else:
         return(False)

def return_sample_names_group(grp):
    """
    given a group, return the names and the column_name associated with that
    """
    compare_dict = load_comparisons()
    for key in compare_dict.keys():
        grp_names = check_key(compare_dict[key],grp)
        if grp_names:
            column_name = compare_dict[key]['column_name'][0]
            return(grp_names,column_name)
    return("","")

def parse_sample_csv_majiq(sample_csv_path):
    """this function takes in the sample csv path and returns a string of
    group and samples with whatever the bam suffix is for the processed bams
    for building a majiq config file from the standard sample csv file used
    by the rna seq pipeline
    """
    samples = pd.read_csv(config['sample_csv_path'])
    #there should be a column which allows you to exclude samples
    samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
    #we're making a dictionary here where each group is the key and all the sample_names
    #in that group are the items then we append the bam_suffix from the config
    temp_dict = (samples2.groupby(config['group_column'])
              .apply(lambda x: set(x['sample_name'] + config['bam_suffix']))
              .to_dict())
    #now we make an empty list, and then fill it up with the values in the dictionary
    #so that the end result looks like group=samplename_bam_suffix,samplename_bam_suffix, etc
    conditions_bams_parsed = []
    for key in temp_dict:
        conditions_bams_parsed.append(key + "=" + ','.join(temp_dict[key]))

    return(conditions_bams_parsed)

def return_parsed_extra_params(extra_params):
    """this function takes extra parameters for a given tool from the config file and parses it
    so that it can be stuck to the end of any shell command
    """
    #starting blank
    if extra_params is None:
        return("")
    cmd = ""
    #for key in extra parameters
    for key in extra_params:
        #append the key value pair if it's a parmeter that needed something
        if extra_params[key]:
            cmd += " --{0} {1}".format(key,extra_params[key])
        else: #otherwise if it's parameter that's just a flag, append just the flag
            cmd += " --{0}".format(key)
    return(cmd)
def majiq_files_by_group(grp):
    """
    given a particular group, return all the majiq files for that group using the samples csv file
    """
    samples = pd.read_csv(config['sample_csv_path'])
    #there should be a column which allows you to exclude samples
    samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
    #read in the comparisons and make a dictionary of comparisons, comparisons needs to be in the config file
    compare_dict = load_comparisons()

    majiq_files = [os.path.join(config['majiq_top_level'],"builder",x) for x in list(samples2.loc[samples2.group == grp].sample_name + config['bam_suffix'] + ".majiq")]

    return(majiq_files)

def majiq_files_from_contrast(grp):
    """
    given a contrast name or list of groups return a list of the files in that group
    """
    #reading in the samples
    samples = pd.read_csv(config['sample_csv_path'])
    #there should be a column which allows you to exclude samples
    samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
    #read in the comparisons and make a dictionary of comparisons, comparisons needs to be in the config file
    compare_dict = load_comparisons()
    #go through the values of the dictionary and break when we find the right groups in that contrast
    grps, comparison_column = return_sample_names_group(grp)
    #take the sample names corresponding to those groups
    if comparison_column == "":
        print(grp)
        return([""])
    grp_samples = list(set(list(samples2[samples2[comparison_column].isin(grps)].sample_name)))

    #build a list with the full path from those sample names
    majiq_files = [os.path.join(config['majiq_top_level'],"builder",x + config['bam_suffix'] + ".majiq") \
                   for x in grp_samples]
    majiq_files = list(set(majiq_files))
    return(majiq_files)

def return_bases_and_contrasts():
    """
    returns all the bases and contrasts from the comparisons.yaml
    """
    comparisons = "config/comparisons.yaml"
    with open(comparisons, 'r') as stream:
        try:
            compare_dict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    base_keys = []
    contrast_keys = []

    for key in compare_dict:
        temp = compare_dict[key]
        for ind, k2 in enumerate(temp):
            if ind == 1:
                base_keys.append(k2)
            if ind == 2:
                contrast_keys.append(k2)
    return(base_keys,contrast_keys)

def get_single_psi_parsed_files():
    """
    return a list of files that will exist
    """
    samples = pd.read_csv(config['sample_csv_path'])
    #there should be a column which allows you to exclude samples
    samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]

    parsed_psi_files = [os.path.join(config['majiq_top_level'],"psi_voila_tsv_single",x) for x in list(samples2.sample_name + config['bam_suffix'] + "_parsed.csv")]

    return(parsed_psi_files)

# takes the featurcounts strand and returns the interpretation for kallisto_output_folder
def get_cuff_strand(fcStrand):
    if fcStrand == "none":
        return("fr-unstranded")
    elif fcStrand == "forward":
        return("fr-secondstrand")
    elif fcStrand == "reverse":
        return("fr-firststrand")


# def get_lib_size(sample_name):
#     """
#     reads the star log files to find the library size for SGSeq
#     """
# getLibSize <- function(logFile){
#       """
#       returns all the bases and contrasts from the comparisons.yaml
#       """
#     stopifnot(file.exists(logFile))
#     log <- readLines(logFile)
#     unique <- log[ grepl("Uniquely mapped reads number", log)]
#     multi <- log[ grepl("Number of reads mapped to multiple loci", log)]
#
#     num_unique <- str_trim( str_split_fixed( unique, "\\|\t", 2)[,2] )
#     num_multi <- str_trim( str_split_fixed( multi, "\\|\t", 2)[,2] )
#     libSize <- as.numeric(num_unique) + as.numeric(num_multi)
#     return(libSize)
#  }
