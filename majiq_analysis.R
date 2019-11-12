library(pacman)
p_load(data.table,janitor)
#read in the majiq tsv files for both the wt vs. f210i and the wt vs wt arsenite treatment 
#for easy matching with stability strip off the last part of the gene name
wt_f210i_majiq = as.data.table(clean_names(fread("/Users/annaleigh/cluster/alb_projects/data/4su_tdp_f210i/majiq2/wt_f21oi_all.tsv")))
wt_f210i_majiq[,gene_simple := tstrsplit(gene_id, "\\.")[1]]
wt_ars = as.data.table(clean_names(fread("/Users/annaleigh/cluster/alb_projects/data/4su_tdp_f210i/majiq2/wtBaseArs_wtDuringArs.tsv")))
wt_ars[,gene_simple := tstrsplit(gene_id, "\\.")[1]]
f210i_ars = as.data.table(clean_names(fread("/Users/annaleigh/cluster/alb_projects/data/4su_tdp_f210i/majiq2/f210iBase_f210iDuringArs.tsv")))
f210i_ars[,gene_simple := tstrsplit(gene_id, "\\.")[1]]


wt_base_post = as.data.table(clean_names(fread("/Users/annaleigh/cluster/alb_projects/data/4su_tdp_f210i/majiq2/wtBaseArse_wtPostArs.tsv")))
wt_base_post[,gene_simple := tstrsplit(gene_id, "\\.")[1]]


wt_during_post = as.data.table(clean_names(fread("/Users/annaleigh/cluster/alb_projects/data/4su_tdp_f210i/majiq2/wtDuringArss_wtPostArse.tsv")))

gsco <- getGScores("phastCons60way.UCSC.mm10")


ir_f210 = wt_f210i_majiq[ir_coords != ""]
ir_f210[,c("start","end") := tstrsplit(ir_coords,"-")]
ir_f210_gr = makeGRangesFromDataFrame(ir_f210)

ir_f210_phylo = gscores(gsco, ir_f210_gr)

ir_wt_ars = wt_ars[ir_coords != ""]
ir_wt_ars[,c("start","end") := tstrsplit(ir_coords,"-")]
ir_wt_ars_gr = makeGRangesFromDataFrame(ir_wt_ars)

wt_ir_phylo = gscores(gsco, ir_wt_ars_gr)






#reading in the data on stability
wide_new_fractions = fread("/Users/annaleigh/Documents/GitHub/4su_tagging/data/wide_bayes_extra_information.csv")
shared_splice_genes = unlist(tstrsplit(intersect(wt_ars$lsv_id, wt_f210i_majiq$lsv_id), "\\.")[1])
gost_shared_splice = gprofiler2::gost(shared_splice_genes, "mmusculus")
gprofiler2::gostplot(gost_shared_splice)
tf_regex = "(?<=Factor: )(.*)(?=; motif)"

gprofiler2::gostplot(gprofiler2::gost(sort(wt_ars[,unique(number_gene_name)]), "mmusculus"))
