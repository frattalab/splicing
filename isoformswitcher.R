
#ensure you have 
if (!require("pacman")) install.packages("pacman")
library(pacman)

p_load("IsoformSwitchAnalyzeR")

top_level = "/Users/annaleigh/cluster/alb_projects/data/muscle/analysis/kallisto"
grep_pattern = "Ctrl[[:digit:]]|SBMA_[[:digit:]]"

baseline = "Ctrl"
comparison = "SBMA"

kallisto_control_sbma = file.path(list.files(top_level,pattern = grep_pattern, full.names = TRUE),"abundance.tsv")
names(kallisto_control_sbma) = list.files(top_level,pattern = grep_pattern) 
k_quant_control_sbma = importIsoformExpression(sampleVector = kallisto_control_sbma,addIsofomIdAsColumn = TRUE)


ctrl_sbma_des <- data.frame(
  sampleID = colnames(k_quant_control_sbma$abundance)[-1],
  condition = ifelse(grepl(baseline,colnames(k_quant_control_sbma$abundance)[-1]),baseline, comparison)
)

exons = "/Users/annaleigh/cluster/vyplab_reference_genomes/annotation/human/ensembl/Homo_sapiens.GRCh38.97.chr_patch_hapl_scaff.gtf.gz"
ntfasta = "/Users/annaleigh/cluster/vyplab_reference_genomes/sequence/human/ensembl/Homo_sapiens.GRCh38.cdna.all.fa.gz"
#exons = "/Users/annaleigh/cluster/vyplab_reference_genomes/annotation/mouse/ensembl/Mus_musculus.GRCm38.97.chr_patch_hapl_scaff.gtf.gz"
#ntfasta = "/Users/annaleigh/cluster/vyplab_reference_genomes/sequence/mouse/ensembl/Mus_musculus.GRCm38.cdna.all.fa.gz"
### Create switchAnalyzeRlist

SwitchList <- importRdata(
  isoformCountMatrix   = k_quant_control_sbma$counts,
  isoformRepExpression = k_quant_control_sbma$abundance,
  designMatrix         = ctrl_sbma_des,
  isoformExonAnnoation = exons,
  isoformNtFasta       = ntfasta,
  removeNonConvensionalChr = TRUE
)


SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = SwitchList,
  geneExpressionCutoff = 10,
  isoformExpressionCutoff = 1,
  removeSingleIsoformGenes = TRUE
)

SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = SwitchListFiltered
)

SwitchListAnalyzed <- analyzeORF(
  SwitchListAnalyzed,
  orfMethod = "longest",
  showProgress=FALSE
)

SwitchListAnalyzed <- extractSequence(
  SwitchListAnalyzed, 
  pathToOutput = '/Users/annaleigh/Documents/data/isoform_switcher/control_sbma/'
)

SwitchListAnalyzed = analyzeCPC2(SwitchListAnalyzed, 
                                    pathToCPC2resultFile = "/Users/annaleigh/Documents/data/isoform_switcher/control_sbma/result_cpc2.txt",
                                    removeNoncodinORFs = FALSE)

SwitchListAnalyzed = analyzePFAM(SwitchListAnalyzed, 
                                    pathToPFAMresultFile = "/Users/annaleigh/Documents/data/isoform_switcher/control_sbma/pfam_annotation.txt")

SwitchListAnalyzed = analyzeSignalP(SwitchListAnalyzed, 
                                       pathToSignalPresultFile = "/Users/annaleigh/Documents/data/isoform_switcher/control_sbma/output_protein_type.txt")


SwitchListAnalyzed <- analyzeAlternativeSplicing(SwitchListAnalyzed, quiet=TRUE)

splicingEnrichment <- extractSplicingEnrichment(
  SwitchListAnalyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)
extractSwitchSummary(SwitchListAnalyzed)

consequencesToAnalyze=c(
  'intron_retention',
  'coding_potential',
  'ORF_seq_similarity',
  'NMD_status',
  'domains_identified')

exampleSwitchListAnalyzed <- analyzeSwitchConsequences(
  SwitchListAnalyzed,
  onlySigIsoforms = TRUE,
  consequencesToAnalyze = consequencesToAnalyze
)




switching_genes = extractTopSwitches(exampleSwitchListAnalyzed, sortByQvals = TRUE, n = 50)$gene_name

pdf("/Users/annaleigh/Documents/data/isoform_switcher/control_sbma/isoform_switches.pdf")
for(g in switching_genes){
  switchPlot(SwitchListAnalyzed, gene = g)
}
dev.off()
