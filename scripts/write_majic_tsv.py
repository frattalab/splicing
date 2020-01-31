#!/usr/bin/python

import argparse,csv, yaml
CONFIGFILE = "/SAN/vyplab/alb_projects/pipelines/splicing/config/config.yaml"
def parse_sample_csv_majiq(sample_csv_path):
    """this function takes in the sample csv path and returns a string of
    group and samples with whatever the bam suffix is for the processed bams
    for building a majiq config file from the standard sample csv file used
    by the rna seq pipeline
    """
    with open(CONFIGFILE, 'r') as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    samples = pd.read_csv(sample_csv_path)
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

def write_tsv_file(tsv_output_destination):
        with open(CONFIGFILE, 'r') as stream:
            try:
                config = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        conditions_bams_parsed = parse_sample_csv_majiq(config['sample_csv_path'])
        options = [
            "[info]",
            "readlen=" + str(config['read_len']),
            "bamdirs=" + config['bam_dir'],
            "genome=" + config['genome_refcode'],
            "strandness=" + config['strand_code'],
            "[experiments]"]

        options += conditions_bams_parsed
        with open(tsv_output_destination,"w") as outfile:
            for opt in options:
                outfile.write(opt + "\n")
        return()


def main(argv):
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                        help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                        const=sum, default=max,
                        help='sum the integers (default: find the max)')

    args = parser.parse_args()
    print(args.accumulate(args.integers))

if __name__ == "__main__":
   main()
