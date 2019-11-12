p_load(IsoformSwitchAnalyzeR)
kallisto_wt_ars = file.path(list.files("/Users/annaleigh/cluster/alb_projects/data/4su_tdp_f210i/kallisto",pattern = "WT[1,2,3][A,B]", full.names = TRUE),"abundance.tsv")
names(kallisto_wt_ars) = list.files("/Users/annaleigh/cluster/alb_projects/data/4su_tdp_f210i/kallisto",pattern = "WT[1,2,3][A,B]") 
k_quant_wt_ars = importIsoformExpression(sampleVector = kallisto_wt_ars,addIsofomIdAsColumn = TRUE)


wt_ars_des <- data.frame(
  sampleID = colnames(k_quant_wt_ars$abundance)[-1],
  condition = ifelse(grepl("B",colnames(k_quant_wt_ars$abundance)[-1]),"during", "base"))


exons = "/Users/annaleigh/cluster/vyplab_reference_genomes/annotation/mouse/ensembl/Mus_musculus.GRCm38.97.chr_patch_hapl_scaff.gtf.gz"
ntfasta = "/Users/annaleigh/cluster/vyplab_reference_genomes/sequence/mouse/ensembl/Mus_musculus.GRCm38.cdna.all.fa.gz"
### Create switchAnalyzeRlist

arsSwitchList <- importRdata(
  isoformCountMatrix   = k_quant_wt_ars$counts,
  isoformRepExpression = k_quant_wt_ars$abundance,
  designMatrix         = wt_ars_des,
  isoformExonAnnoation = exons,
  isoformNtFasta       = ntfasta,
  removeNonConvensionalChr = TRUE
)



arsSwitchListFiltered <- preFilter(
  switchAnalyzeRlist = arsSwitchList,
  geneExpressionCutoff = 10,
  isoformExpressionCutoff = 3,
  removeSingleIsoformGenes = TRUE
)

arsSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = arsSwitchListFiltered
)

arsSwitchListAnalyzed <- analyzeORF(
  arsSwitchListAnalyzed,
  orfMethod = "longest",
  # genomeObject = Hsapiens, # not necessary since sequences are already stored in the switchAnalyzeRlist
  showProgress=FALSE
)

arsSwitchListAnalyzed <- extractSequence(
  arsSwitchListAnalyzed, 
  pathToOutput = '/Users/annaleigh/Documents/data/isoform_switcher/wt_ars/'
)

arsSwitchListAnalyzed = analyzeCPC2(arsSwitchListAnalyzed, 
                                    pathToCPC2resultFile = "/Users/annaleigh/Documents/data/isoform_switcher/wt_ars/result_cpc2.txt",
                                    removeNoncodinORFs = FALSE)
arsSwitchListAnalyzed = analyzePFAM(arsSwitchListAnalyzed, 
                                    pathToPFAMresultFile = "/Users/annaleigh/Documents/data/isoform_switcher/wt_ars/pfam_annotation.txt")

arsSwitchListAnalyzed = analyzeSignalP(arsSwitchListAnalyzed, 
                                       pathToSignalPresultFile = "/Users/annaleigh/Documents/data/isoform_switcher/wt_ars/output_protein_type.txt")


arsSwitchListAnalyzed <- analyzeAlternativeSplicing(arsSwitchListAnalyzed, quiet=TRUE)

splicingEnrichment <- extractSplicingEnrichment(
  arsSwitchListAnalyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)
extractSwitchSummary(arsSwitchListAnalyzed)

consequencesToAnalyze=c(
  'intron_retention',
  'coding_potential',
  'ORF_seq_similarity',
  'NMD_status',
  'domains_identified')

arsSwitchListAnalyzed <- analyzeSwitchConsequences(
  arsSwitchListAnalyzed,
  onlySigIsoforms = TRUE,
  consequencesToAnalyze = consequencesToAnalyze
)




switching_genes = extractTopSwitches(arsSwitchListAnalyzed, filterForConsequences = TRUE, sortByQvals = TRUE, n = 50)$gene_name

pdf("/Users/annaleigh/Documents/data/isoform_switcher/wt_ars/swithces.pdf")
for(g in switching_genes){
  switchPlot(arsSwitchListAnalyzed, gene = g)
}
dev.off()
library(gridExtra)

for (i in seq(length(plot_list))) {
  do.call("grid.arrange", p[[i]])  
}
dev.off()
