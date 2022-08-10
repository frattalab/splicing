
library(GenomicFeatures)
# gtf_path = '/SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v38.annotation.gtf'
# gtf_path = "/Users/annaleigh/cluster/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v38.annotation.gtf"
gtf_path = '/Users/annaleigh/Downloads/gencode.v38.annotation.gtf.gz'
print(stringr::str_c(Sys.time(), " - Importing the GTF..."))
gtf_obj =  GenomicFeatures::makeTxDbFromGFF(gtf_path,format = 'gtf')
g = as.data.table(rtracklayer::import(gtf_path))
biotype_table = unique(g[,.(gene_id,gene_type,transcript_id,transcript_type,gene_name)])
biotype_table = biotype_table %>% 
    mutate(gene_id = gsub("\\..*","",gene_id)) %>%
    mutate(transcript_id = gsub("\\..*","",transcript_id))
rm(g)
beepr::beep(6)

top_folder = '/Users/annaleigh/Documents/GitHub/tdp_43_psi_rankings/'
experiment = 'chx/CycloheximideControl-CycloheximideTDP43KD'


transcript_counts = fread(glue::glue("{top_folder}transcript_counts/{experiment}.transcript_counts.csv"))
transcript_counts = transcript_counts %>% 
    mutate(gene_id = gsub("\\..*","",GENEID)) %>%
    mutate(transcript_id = gsub("\\..*","",TXNAME))

casettes = rtracklayer::import(glue::glue('{top_folder}modules_output/{experiment}.cryptic_exons.bed'))

#get regions per transcript

tx = transcripts(gtf_obj)

#each exon now which transcript it could fall into
annotated_casettes = as.data.table(annotatr::annotate_regions(sort(casettes),tx))
annotated_casettes = annotated_casettes[strand == annot.strand]
annotated_casettes = annotated_casettes %>% 
    mutate(transcript_id = gsub("\\..*","",annot.tx_name))

#find the overall transcript expression in this dataset

transcript_expression = transcript_counts %>% 
    melt() %>% 
    group_by(transcript_id) %>% 
    summarize(avg_expression = mean(value),gene_id) %>% 
    ungroup() %>% 
    unique() %>% 
    mutate()

protein_coding_genes = annotated_casettes %>% 
    left_join(transcript_expression,by = 'transcript_id') %>% 
    filter(!is.na(avg_expression)) %>%  #get rid of transcripts with no expression 
    unique() %>% 
    left_join(biotype_table, by = c('transcript_id','gene_id')) %>% 
    filter(gene_type == "protein_coding") %>% 
    filter(transcript_type == "protein_coding")

max_tx_protein_coding_genes = protein_coding_genes%>% 
    mutate(cryptic_coords = paste0(seqnames,":",start,"-",end)) %>% 
    group_by(cryptic_coords,gene_id) %>% 
    mutate(max_by_gene = max(avg_expression)) %>% 
    ungroup() %>% 
    filter(avg_expression == max_by_gene) %>% 
    dplyr::select(seqnames:strand,avg_expression,gene_name,gene_id,transcript_id,cryptic_coords) %>% 
    dplyr::select(-width) %>% 
    unique()


fwrite(max_tx_protein_coding_genes, glue::glue('{top_folder}modules_output/{experiment}.cryptic_exons.protein_coding_with_transcript_backbone.csv'))
