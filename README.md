# splicing
Splicing done with MAJIQ tool - This is very bare bones implentation as it stands is a **work in progress**, 

The purpose of this pipeline is to be able to run MAJIQ using Snakemake. The aim is to make MAJIQ easier to run for non-bioinformaticians and it produces additional parsing and annotation to the MAJIQ output. 

# Needed files
1. Aligned, sorted, and indexed BAM files of RNA-seq. You will need .bam and .bai files for all your samples. 
2. GFF3 and GTF of your species of interest
3. A formatted sample sheet, see examples and explanation below
# Final outputs
Underneath the folder in 

`majiq_top_level: /SAN/vyplab/alb_projects/data/linked_bams_f210i_brain/majiq/`
```
majiq
├── builder
│   ├── wt_sample1.majiq
│   ├── wt_sample1.sj
│   ├── wt_sample2.majiq
│   ├── wt_sample2.sj
│   ├── mut_sample1.majiq
│   ├── mut_sample1.sj
│   ├── mut_sample2.majiq
│   ├── mut_sample2.sj
│   ├── majiq.log
│   └── splicegraph.sql
├── delta_psi
│   ├── wt_mut.deltapsi.tsv
│   ├── wt_mut.deltapsi.voila
│   └── deltapsi_majiq.log
├── delta_psi_voila_tsv
│   ├── wt_mut.junctions.bed
│   ├── wt_mut.csv
│   ├── wt_mut.gff3
│   ├── wt_mut_parsed_psi.tsv
│   └── wt_mut.psi.tsv
├── run_name_majiqConfig.tsv
├── psi_single
│   ├── wt_sample1.tsv
│   ├── wt_sample1.voila
│   ├── wt_sample2.tsv
│   ├── wt_sample2.voila
│   ├── mut_sample1.tsv
│   ├── mut_sample1.voila
│   ├── mut_sample2.tsv
│   ├── mut_sample2.voila
├── psi_voila_tsv_singlehis
│   ├── wt_sample1.tsv
│   ├── wt_sample1.voila
│   ├── wt_sample2.tsv
│   ├── wt_sample2.voila
│   ├── mut_sample1.tsv
│   ├── mut_sample1.voila
│   ├── mut_sample2.tsv
│   ├── mut_sample2.voila
└── psi
    ├── wt.psi.tsv
    ├── wt.psi.voila
    ├── mut.psi.tsv
    ├── mut.psi.voila
    └── psi_majiq.log
```
## Submitting on SGE

1. Build step
`source submit.sh build run_name`
2. PSI step
`source submit.sh psi run_name`

with whatever run name you'd like

## Submitting on Slurm

1. Build step
`source submit_slurm.sh build run_name`
2. PSI step
`source submit_slurm.sh psi run_name`
with whatever run name you'd like
## Running without a cluster

If you don't have a cluster, you can run straight with snakemake
`snakemake -s workflows/build.smk`
`snakemake -s workflows/psi.smk`



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

Please cite Dasper, Snakemake, and of course MAJIQ if you use this pipeline.
