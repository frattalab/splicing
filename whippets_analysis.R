wt_f_base = fread("/Users/annaleigh/Documents/data/whippet_snake/wt_f.diff.gz")

corrected_names = names(wt_f_base)[2:length(wt_f_base)]
wt_f_base[,Entropy := NULL]
names(wt_f_base) = corrected_names
wide_new_fractions = fread("/Users/annaleigh/Documents/GitHub/4su_tagging/data/wide_bayes_extra_information.csv")
wide_new_fractions_base = wide_new_fractions[cond == "base"]
setDT(wt_f_base)[wide_new_fractions_base, bayesian_p := bayesian_p, on = c(Gene = "gene")]
setDT(wt_f_base)[wide_new_fractions_base, median_diff := median_diff, on = c(Gene = "gene")]
setDT(wt_f_base)[wide_new_fractions_base, mean_s1_ntr := mean_s1_ntr, on = c(Gene = "gene")]
setDT(wt_f_base)[wide_new_fractions_base, mean_f210_ntr := mean_s2_ntr, on = c(Gene = "gene")]
#exon skipping cassette exon usage even is CE 
probability = 0.95 # min probability
psiDelta = 0.1 # min change in PSI
eventTypes = "all" # all event types
minCounts = 100 # mean of at least 100 counts in one condition
wt_mut_base = readWhippetDataSet("/Users/annaleigh/Documents/data/whippet_snake/")

whippet_sampleTable <- data.frame(sample=c("F2A_1","F3A_1","WT1A_1","WT2A_1","WT3A_1"),
                                  condition=c("f","f","wt","wt","wt"),
                                  replicate=(c("1","2","1","2","3")))

wds_filters <- filterWhippetEvents(
  whippetDataSet = wt_mut_base,
  probability = 0.95, # min probability
  psiDelta = 0.1, # min change in PSI
  eventTypes = "all", # all event types
  minCounts = 100, # mean of at least 100 counts in one condition
  sampleTable = whippet_sampleTable)

dif_spli = as.data.table(wds_filters@diffSplicingResults)
