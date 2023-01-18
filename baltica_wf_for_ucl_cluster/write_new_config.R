library(tidyverse)
library(yaml)

config <- dplyr::lst(
    # path to the project directory, where the output files
    # are written
    path = "/home/tbrittoborges/sirv_test",

    # path to where the data files are
    sample_path = "/home/tbrittoborges/Baltica/data",

    # assembly identifier is critical to generate the genome browser links
    assembly = "SIRV", # this is artificial, so no working genome browser links

    # path to the reference transcritome annotation, gtf format
    ref = "/home/tbrittoborges/Baltica/data/SIRV.gtf",

    # path to the reference transcritome sequence, fasta format
    ref_fa = "/home/tbrittoborges/Baltica/data/SIRV.gtf",

    # What library preparation did you use?
    # see https://chipster.csc.fi/manual/library-type-summary.html
    # dUTP: reverse
    # ligation: foward
    # unstranded: unstranded
    strandness = "reverse",

    # maximum read lenght after trimming
    read_len = 150L,

    # method specific parameters

    #
    leafcutter_min_samples_per_group = 2L,
    leafcutter_min_samples_per_intron = 2L
)


files <- Sys.glob(file.path(config$sample_path, "*.bam"))
# checking if we get all the files
length(files)

# here we create groups for each file,
# specifically using the file names,
# but this can also be done manually
groups <- case_when(
    grepl("mix1", files) ~ "mix1",
    grepl("mix2", files) ~ "mix2",
    grepl("mix3", files) ~ "mix3"
)
groups

# baltica uses the {group}_{rep} format
# to identify replicates of the same experiment
groups_reps <- paste(
    groups,
    ave(rep(NA, length(groups)), groups, FUN = seq_along),
    sep = "_"
)
groups_reps

# What comparisons one want to make?
# here we compare the three groups
# all versus all
contrasts <- list(
    c("mix2", "mix1"),
    c("mix3", "mix1"),
    c("mix3", "mix2")
)
names(contrasts) <- sapply(contrasts, paste0, collapse = "-vs-")

samples <- as.list(setNames(nm = groups_reps, files))
head(samples)

config[["contrasts"]] <- contrasts
config[["samples"]] <- samples

# check the config before writting
config

# by default, we write the config where the data is
yaml::write_yaml(config, file.path(config$sample_path, "config.yaml"))

yaml::write_yaml(
    list(
        "__default__" = list(
            "mem" = 16000L,
            "cpu" = 4L,
            "out" = "logs/{rule}_{wildcards}.log",
            "partition" = "general"
        ),
        "rseqc_read_duplication" = list(
            "mem" = 64000L
        ),
        "build" = list(
            "mem" = 6400L,
            "cpu" = 4L,
            "partition" = "long"
        ),
        "qc" = list(
            "mem" = 32000L
        ),
        "junctioseq_analysis" = list(
            "mem" = 64000L,
            "cpu" = 10L
        ),
        "run_rmats" = list(
            "mem" = 64000L
        )
    ),
    file.path(base_path, "cluster.yaml")
)