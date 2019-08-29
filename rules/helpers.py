import os
import pandas as pd

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
    cmd = ""
    #for key in extra parameters
    for key in extra_params:
        #append the key value pair if it's a parmeter that needed something
        if extra_params[key]:
            cmd += " --{0} {1}".format(key,extra_params[key])
        else: #otherwise if it's parameter that's just a flag, append just the flag
            cmd += " --{0}".format(key)
    return(cmd)

def majiq_files_by_group(group):
    """
    given a particular group, return all the majiq files for that group using the samples csv file
    """
    samples = pd.read_csv(config['sample_csv_path'])
    #there should be a column which allows you to exclude samples
    samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
    majiq_files = [os.path.join(config['majiq_top_level'],"builder",x) for x in list(samples2.loc[samples2.group == group].sample_name + config['bam_suffix'] + ".majiq")]
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
    return([base_keys,contrast_keys])
