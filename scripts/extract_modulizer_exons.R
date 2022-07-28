library(data.table)
library(tidyverse)
library(rlang)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringi)

gtf_path = "/Users/annaleigh/cluster/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v38.annotation.gtf"
print(stringr::str_c(Sys.time(), " - Creating Exon by Transcript DT..."))
gtf_obj =  GenomicFeatures::makeTxDbFromGFF(gtf_path,format = 'gtf')



# now read in the annotated GTF from gencode
gtf_exons = as.data.table(GenomicFeatures::exons(gtf_obj,use.names = TRUE))
gtf_exons = gtf_exons %>% mutate(coords = paste0(seqnames,":",start,"-",end))
    
exon_by_transcript = unlist(exonsBy(gtf_obj,'tx',use.names = TRUE))









process_event_files = function(path_to_file, event_name = 'alt3_5'){
    cryps = cryptic_junctions %>% 
                 pull(paste_into_igv_junction)
    
    the_event  = fread(path_to_file)
    the_event[,paste_into_igv_junction := paste0(seqid,":",junction_coord)]
    the_event[,spliced_with_exon := paste0(seqid,":",spliced_with_coord)]
    the_event[,reference_exon := paste0(seqid,":",reference_exon_coord)]
    
    cryptic_exon_of_event = the_event %>% 
        filter(paste_into_igv_junction %in% cryps) %>% 
        select(gene_name,contains("psi"),spliced_with_exon,reference_exon,paste_into_igv_junction,strand) %>% 
        mutate(type = event_name) %>% 
        unique()
    
    return(cryptic_exon_of_event)
}

top_folder = '/Users/annaleigh/cluster/alb_projects/data/tdp_ko_collection/invitro_majiq/majiq/'
bottom_folder = 'controlSeddighiCorticalNeuron-TDP43KDSeddighiCorticalNeuron'

#top_folder = "/Users/annaleigh/cluster/first_weeks/TDP_CHX_CLONES_GLIA/chx/majiq/"
#bottom_folder = 'CycloheximideControl-CycloheximideTDP43KD'

baseline = str_split(bottom_folder,"-")[[1]][1]
contrast = str_split(bottom_folder,"-")[[1]][2]
basepsi = paste0(baseline, '_median_psi')
contrastpsi = paste0(contrast, '_median_psi')

prob_column = paste0(bottom_folder,"_probability_changing")
# read in the annotated junctions -----------------------------------------
middle_folder = 'delta_psi_voila_tsv' 
my_junctions = fread(file.path(top_folder,middle_folder,paste0(bottom_folder,'_annotated_junctions.csv')))
my_junctions = my_junctions %>% 
    select(seqnames,gene_id,gene_name,mean_dpsi_per_lsv_junction:lsv_type,
           paste_into_igv_junction,junc_in_ref,junc_cat,strand) %>% 
    select(-lsv_type) 

middle_folder = 'modulizers/' 
junctions_mod = fread(file.path(top_folder,middle_folder,paste0(bottom_folder,'/junctions.tsv')))
junctions_mod = junctions_mod[junction_coord != '']
junctions_mod[,paste_into_igv_junction := paste0(seqid,":",junction_coord)]
junctions_mod = junctions_mod %>% 
    separate(junction_coord, into = c("start","end"),sep = "-",remove = FALSE)

cryptic_junctions = junctions_mod %>% 
    filter(junctions_mod[,18] < 0.05 & junctions_mod[,19] > 0.1) %>% 
    select(4,13,14,5,2,3,6,9,18,19,20,21,23) %>% 
    unique() %>% 
    mutate(module_event_combination = gsub("\\^\\d+","",module_event_combination)) %>% 
    separate_rows(module_event_combination,sep = "\\|") %>% 
    unique() %>% 
    left_join(my_junctions[,.(paste_into_igv_junction,junc_cat)]) %>% 
    unique()

folder_path = file.path(top_folder,middle_folder,paste0(bottom_folder))

read_me = list.files(folder_path)
read_me = read_me[!(read_me %in% c('junctions.tsv','voila.log','heatmap.tsv','orphan_junction.tsv'))]

for(r in read_me){
    print(r)
    
}
print(unique(cryptic_junctions$module_event_combination))

# process the exons -----------------------------------------------------------
alt3_5_exons = process_event_files(file.path(folder_path,'alt3and5prime.tsv'),event_name = 'alt3_5')
alt3_exons = process_event_files(file.path(folder_path,'alt3prime.tsv'),event_name = 'alt3')
alt5_exons = process_event_files(file.path(folder_path,'alt5prime.tsv'),event_name = 'alt5')
afe_exons = process_event_files(file.path(folder_path,'alternate_first_exon.tsv'),event_name = 'afe')
ale_exons = process_event_files(file.path(folder_path,'alternate_last_exon.tsv'),event_name = 'ale')

cassette_exons = process_event_files(file.path(folder_path,'cassette.tsv'),event_name = 'cassette')

mutually_exclusive_exons = process_event_files(file.path(folder_path,'mutually_exclusive.tsv'),event_name = 'mxe')
putative_alt3_exons = process_event_files(file.path(folder_path,'p_alt3prime.tsv'),event_name = 'putative_alt3')
putative_alt5_exons = process_event_files(file.path(folder_path,'p_alt5prime.tsv'),event_name = 'putative_alt5')
putative_afe_exons = process_event_files(file.path(folder_path,'p_alternate_first_exon.tsv'),event_name = 'putative_afe')
putative_ale_exons = process_event_files(file.path(folder_path,'p_alternate_last_exon.tsv'),event_name = 'putative_ale')
# all the exons -----------------------------------------------------------

all_exons = rbind(afe_exons, ale_exons, alt3_5_exons, alt3_exons, alt5_exons,
                  cassette_exons,mutually_exclusive_exons,
                  putative_afe_exons, putative_ale_exons,
                  putative_alt3_exons, putative_alt5_exons) %>% 
    left_join(my_junctions[,.(paste_into_igv_junction,junc_cat)]) %>% 
    filter(junc_cat != 'novel_exon_skip') 

melted_exons = all_exons %>% 
    filter(all_exons[,2] < 0.05) %>% 
    filter(junc_cat != "annotated") %>% 
    select(gene_name,spliced_with_exon,reference_exon,strand,paste_into_igv_junction,type,junc_cat) %>% 
    melt(id.vars = c('paste_into_igv_junction','gene_name','strand','type','junc_cat')) %>% 
    unique()


melted_exons = melted_exons[!(value %in% gtf_exons$coords)]

melted_exons %>% 
    separate(value,
             into = c("seqnames",'res'),sep = ":") %>% 
    mutate(res = gsub("--1","-1",res)) %>% 
    mutate(res = gsub("^-1","1",res)) %>% 
    separate(res, into = c("start","end"),sep = "-",remove = FALSE) %>% 
    mutate(start = as.numeric(start),
           end = as.numeric(end)) %>% 
    filter(is.na(start) | is.na(end))

beepr::beep()
    
# 
#     mutate(name = paste0(gene_name,"|",paste_into_igv_junction)) %>% 
#     separate(value,
#              into = c("seqnames",'start','end')) %>% 
#     filter(start != 1) %>% 
#     filter(end != 1) %>% 
#     unique() %>% 
#     filter(!is.na(as.numeric(end))) %>% 
#     filter(!is.na(as.numeric(start))) %>% 
#     makeGRangesFromDataFrame() %>% 
#     rtracklayer::export('cryptic_exons_again_test.bed')
# 
#     cryptic_junctions