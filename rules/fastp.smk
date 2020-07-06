
configfile: "config/config.yaml"
include: "helpers.py"

fastp_trimmed_output_folder = os.path.join(config["top_level_project_folder"],config["fastp_trimmed_output_folder"])

#make sure the output folder for fastptrimming exists before running anything
os.system("mkdir -p {0}".format())
#read in a samples table
SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)
#so I want this rule to be run ONCE for every fast1, so the wild cards I'm giving are the 'name' of the fastq of the first read
#name = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in SAMPLES['fast1'].tolist()]
name = [strpd.rpartition('/')[2].split(".")[0] for strpd in SAMPLES['fast1'].tolist()]


UNITS = SAMPLES['unit'].tolist()
SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

rule merge_all_trimmed:
    input:
        expand(config["fastp_trimmed_output_folder"] + "{unit}_{name}_R1_trimmed.fastq.gz",zip, unit = UNITS,name = SAMPLE_NAMES),
        expand(config["fastp_trimmed_output_folder"] + "{unit}_{name}_R2_trimmed.fastq.gz" if config["end_type"] == "pe" else [],zip, unit = UNITS,name = SAMPLE_NAMES),
        expand(config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz", name = SAMPLE_NAMES),
        expand(config["merged_fastq_folder"] + "{name}_2.merged.fastq.gz" if config["end_type"] == "pe" else [],name = SAMPLE_NAMES)
    wildcard_constraints:
        name="|".join(SAMPLE_NAMES),
        unit="|".join(UNITS)

if config['end_type'] == "pe":
    rule fastp_trimming:
        input:
        #get the value in the fast1 column
            fastq_file = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = True),
            fastq_file2 = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = False)
        wildcard_constraints:
            name="|".join(SAMPLE_NAMES),
            unit="|".join(UNITS)
        output:
            out_fastqc = config["fastp_trimmed_output_folder"] + "{unit}_{name}_R1_trimmed.fastq.gz",
            out_fastqc2 = config["fastp_trimmed_output_folder"] + "{unit}_{name}_R2_trimmed.fastq.gz",
            fastpjson = config["fastp_trimmed_output_folder"] + "{unit}_{name}_fastp.json",
            fastphtml = config["fastp_trimmed_output_folder"] + "{unit}_{name}_fastp.html",
        params:
            fastp_parameters = return_parsed_extra_params(config['fastp_parameters']),
            fastpjson = config["fastp_trimmed_output_folder"] + "{unit}_{name}_fastp.json",
            fastphtml = config["fastp_trimmed_output_folder"] + "{unit}_{name}_fastp.html"
        shell:
            """
            {config[fastp_path]} --in1 {input.fastq_file} --in2 {input.fastq_file2} --out1 {output.out_fastqc} --out2 {output.out_fastqc2} --json {output.fastpjson} --html {output.fastphtml} {params.fastp_parameters}
            """
else:
        rule fastp_trimming:
            input:
            #get the value in the fast1 column
                fastq_file = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = True)
            wildcard_constraints:
                name="|".join(SAMPLE_NAMES),
                unit="|".join(UNITS)
            output:
                out_fastqc = config["fastp_trimmed_output_folder"] + "{unit}_{name}_R1_trimmed.fastq.gz",
                fastpjson = config["fastp_trimmed_output_folder"] + "{unit}_{name}_fastp.json",
                fastphtml = config["fastp_trimmed_output_folder"] + "{unit}_{name}_fastp.html",
            params:
                fastp_parameters = return_parsed_extra_params(config['fastp_parameters']),
                fastq_file2 = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = False),
            #out_fastqc2 = lambda wildcards: return_fastq2_name(wildcards.name,wildcards.unit),
                fastpjson = config["fastp_trimmed_output_folder"] + "{unit}_{name}_fastp.json",
                fastphtml = config["fastp_trimmed_output_folder"] + "{unit}_{name}_fastp.html"
            shell:
                """
                {config[fastp_path]} -i {input.fastq_file} -o {output.out_fastqc} --json {output.fastpjson} --html {output.fastphtml} {params.fastp_parameters}
                """

if config['end_type'] == "pe":
    rule merge_trimmed:
        input:
            one = lambda wildcards: get_trimmed(wildcards.name)[0],
            two = lambda wildcards: get_trimmed(wildcards.name)[1]
        wildcard_constraints:
            name="|".join(SAMPLE_NAMES)
        output:
            out_one = config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz",
            out_two = config['merged_fastq_folder'] + "{name}_2.merged.fastq.gz"
        params:
            #taking the input files and putting them into a comma separated list
            one = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[0]),
            two = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[1])
        shell:
            """
            mv {params.two} {output.out_two}
            mv {params.one} {output.out_one}
            """
else:
    rule merge_trimmed:
        input:
            one = lambda wildcards: get_trimmed(wildcards.name)[0]
        wildcard_constraints:
            name="|".join(SAMPLE_NAMES)
        output:
            out_one = config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz"
        params:
            #taking the input files and putting them into a comma separated list
            one = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[0]),
        shell:
            """
            mv {params.one} > {output.out_one}
            """
