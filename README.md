# splicing
Splicing done with MAJIQ tool - This is very bare bones implentation as it stands is a **work in progress**, 



1. Build step
source submit.sh build run_name
2. PSI step
source submit.sh psi run_name
with whatever run name you'd like

Change line 26 in in submit_majiq.sh to rules/majiq_psi.smk

And then run submit_majiq.sh

And finally
submit_majiq_parse.sh
 
To run the majiq PSI with each sample as it's own group "submit_majiq_single.sh runname"

The scripts folder as of 2020-02-27 contains a series of half-functional things. Be cautious.

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
