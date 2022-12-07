library(tidyverse)
library(gt)


descrip <- tribble(
  ~method, ~parameter, ~description,
  "rMATS", "annotation",  "Uses annotation and library reads to create a database of splice graphs (+-)",
  "JunctionSeq", "annotation", " Uses annotation to name locus (-)",
  "Leafcutter", "annotation", "Uses annotation and library to selected testable features (introns) (+)",
  "Majiq", "annotation", "Uses annotation and library to selected testable events (+)",
  # 
  "rMATS", "coverage_filter", "None",
  "JunctionSeq", "coverage_filter", "--minCount at `qorts mergeNovelSplices` default: 9",
  "Leafcutter", "coverage_filter", "-m at `leafcutter_cluster_regtools.py`` default: 50",
  "Leafcutter", "coverage_filter", "--min_coverage at `leafcutter_ds.R` default: 20",
  "Majiq", "coverage", "--minreads `majiq build` default: 10)",
  # 
  "rMATS", "replicate_filtering", "None",
  "JunctionSeq", "replicate_filtering", "None",
  "Leafcutter", "replicate_filtering", "--min_samples_per_intron at `leafcutter_ds.R` default: 5)",
  "Leafcutter", "replicate_filtering", "--min_samples_per_group at `leafcutter_ds.R` default: 3",
  "Majiq", "replicate_filtering", "--min-experiments at `majiq build` default: 0.5)",
  # 
  "rMATS", "experimental_design", "Case vs control",
  "JunctionSeq", "experimental_design", "Case vs control, supports covariates",
  "Leafcutter", "experimental_design", "Supports design matrix",
  "Majiq", "experimental_design", "Case control, paired-design"
)

