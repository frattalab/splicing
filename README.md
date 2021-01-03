# splicing
Splicing done with MAJIQ tool - This is very bare bones implentation as it stands is a **work in progress**, 


## Submitting on SGE

1. Build step
`source submit.sh build run_name`
2. PSI step
`source submit.sh psi run_name`
with whatever run name you'd like

## Submitting on Slurm

1. Build step
source submit_slurm.sh build run_name
2. PSI step
source submit_slurm.sh psi run_name
with whatever run name you'd like


The scripts folder as of 2020-11-05 contains semi-functional annotation scripts. Be cautious.

See example data for the formating of sample sheets.
The following columns are mandatory:
sample_name,
unit,
fast1,
fast2,
group
exclude_sample_downstream_analysis

unit, fast1, fast2 can be placeholders, they are to maintain the same sample sheet structure across my RNAseq alignment pipeline and the splicing pipeline

exclude_sample_downstream_analysis should be present, if you want to exclude a sample it should be a 1

## Annotation of splicing events
Annotation is done with a function grabbed directly from source code here:
https://github.com/dzhang32/dasper/

Please cite ^ if you use this pipeline.
