library(data.table)
library(tidyverse)
library(rlang)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringi)



mean_exon_psi <- function(my_df, columnA,columnB){
    

    my_df %>% 
        group_by(paste_into_igv_exon) %>% 
        summarise(gene_name.x,
                  strand.x,
                  base_mean = mean(!! ensym(columnA)),
                  contrast_mean = mean(!! ensym(columnB)), .groups = 'drop' ) %>% 
        dplyr::rename(gene_name = gene_name.x, 
                      strand = strand.x)
    
    
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

cryptic_junctions %>% 
    dplyr::rename(base_psi = 9,contrast_psi = 10,seqnames = seqid) %>% 
    mutate(name = paste0(gene_name,"|",base_psi,"|",contrast_psi)) %>% 
    select(seqnames,start,end,strand,name) %>% 
    unique() %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
    rtracklayer::export('cryptic.junctions.bed')

# Read in the summary to find out what kind of files I'll need for each type -----------------------------------------
middle_folder = 'modulizers/' 

# cassette exons -----------------------------------------------------------
folder_path  = file.path(top_folder,middle_folder, bottom_folder)

cassette = fread(file.path(folder_path,"cassette.tsv"))
cassette[,paste_into_igv_junction := paste0(seqid,":",junction_coord)]
cassette[,paste_into_igv_exon := paste0(seqid,":",spliced_with_coord)]

cryptic_casette <- cryptic_junctions %>% 
    select(-c(9,10)) %>% 
    left_join(cassette[denovo == TRUE], by = 'paste_into_igv_junction') %>% 
    filter(!is.na(paste_into_igv_exon)) %>% 
    mean_exon_psi(.,!!basepsi,!!contrastpsi) %>% 
    mutate(name = paste0(gene_name,"|",base_mean,"|",contrast_mean)) %>% 
    separate(paste_into_igv_exon,into =  c("chr","start","end")) %>% 
    mutate(end = as.numeric(end),
           start = as.numeric(start)) %>% 
    filter((end >= start - 1)) %>% 
    unique() %>% 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE) %>% 
    plyranges::filter(base_mean < 0.05) %>% 
    filter(start != 1)


cryptic_casette %>% rtracklayer::export(paste0(bottom_folder,'_cryptic_casette.bed'))

cryptic_casette$seq<- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                                      cryptic_casette))

cryptic_casette$seqwid<- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                                              cryptic_casette + 2))

cryptic_casette$acceptor_motif = substr(cryptic_casette$seqwid,1,2)
cryptic_casette$donor_motif = stri_sub(cryptic_casette$seqwid,-2,-1)

cryptic_casette$seqwid = NULL

cryptic_casette %>% 
    plyranges::mutate(donor_canonical = (donor_motif == "GT"))  %>% 
    plyranges::mutate(acceptor_canonical = (acceptor_motif == "AG"))  %>% 
    plyranges::filter(!donor_canonical)

# last exons -----------------------------------------------------------
last_exons = fread(file.path(folder_path,"alternate_last_exon.tsv"))
last_exons[,paste_into_igv_junction := paste0(seqid,":",junction_coord)]
last_exons[,paste_into_igv_exon := paste0(seqid,":",spliced_with_coord)]


cryptic_exon_junctions %>% 
    left_join(last_exons, by = 'paste_into_igv_junction') %>% 
    filter(!is.na(paste_into_igv_exon))  %>% 
    select(!!basepsi)
    mean_exon_psi(.,!!basepsi,!!contrastpsi) %>% 
    mutate(name = paste0(gene_name,"|",base_mean,"|",contrast_mean)) %>% 
    unique() %>% 
    separate(paste_into_igv_exon,into =  c("chr","start","end")) %>% 
    mutate(end = as.numeric(end),
           start = as.numeric(start)) %>% 
    filter((end >= start - 1)) %>% 
    unique() %>% 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE) 

# putative last exons -----------------------------------------------------------
p_last_exons = fread(file.path(folder_path,"p_alternate_last_exon.tsv"))
p_last_exons[,paste_into_igv_junction := paste0(seqid,":",junction_coord)]
p_last_exons[,paste_into_igv_exon := paste0(seqid,":",spliced_with_coord)]


p_cryptic_last_exons = cryptic_exon_junctions %>% 
    left_join(p_last_exons, by = 'paste_into_igv_junction') %>% 
    filter(!is.na(paste_into_igv_exon))  %>% 
    mean_exon_psi(.,!!basepsi,!!contrastpsi) %>% 
    mutate(name = paste0(gene_name,"|",base_mean,"|",contrast_mean)) %>% 
    unique() %>% 
    separate(paste_into_igv_exon,into =  c("chr","start","end")) %>% 
    mutate(end = as.numeric(end),
           start = as.numeric(start)) %>% 
    filter((end >= start - 1)) %>% 
    unique() %>% 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE) 
