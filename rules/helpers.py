import os
import pandas as pd

def parse_sample_csv_majiq(sample_csv_path):
    """this function takes in the sample csv path and returns a string of
    group and samples with whatever the bam suffix is for the processed bams
    for building a majiq config file from the standard sample csv file used
    by the rna seq pipeline
    """
    samples = pd.read_csv(sample_csv_path)
    #there should be a column which allows you to exclude samples
    samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
    #we're making a dictionary here where each group is the key and all the sample_names
    #in that group are the items then we append the bam_suffix from the config
    temp_dict = (samples2.groupby('group')
              .apply(lambda x: set(x['sample_name'] + config['bam_suffix']))
              .to_dict())
    #now we make an empty list, and then fill it up with the values in the dictionary
    #so that the end result looks like group=samplename_bam_suffix,samplename_bam_suffix, etc
    conditions_bams_parsed = []
    for key in test:
        conditions_bams_parsed.append(key + "=" + ','.join(test[key]))

    return(conditions_bams_parsed)
