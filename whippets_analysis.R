wt_f_base = fread("/Users/annaleigh/Documents/data/whippet_snake/wt_f_baseline.diff.gz")

corrected_names = names(wt_f_base)[2:length(wt_f_base)]
wt_f_base[,Entropy := NULL]
names(wt_f_base) = corrected_names

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
