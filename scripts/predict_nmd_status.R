library(stringi)
library(msa)
library(plyranges)
my_orf = function(transcripts, BSgenome = NULL, returnLongestOnly = TRUE, 
                  allFrames = FALSE,  
                  uORFs = FALSE) 
{
    if (allFrames == TRUE) {
        returnLongestOnly = FALSE
        longest = 1
    }
    browser()
    transcripts$exon_number <- as.numeric(transcripts$exon_number)
    order <- order(transcripts$transcript_id, transcripts$exon_number)
    transcripts <- transcripts[order]
    transcripts$seq <- as.character(Biostrings::getSeq(BSgenome, 
                                                       transcripts))
    seqCat <- aggregate(seq ~ transcript_id, mcols(transcripts), 
                        function(x) (paste(x, collapse = "")))
    ids <- as.character(seqCat$transcript_id)
    
    seqCat <- seqCat$seq
    
    rm <- which(grepl("N", seqCat))
    if (length(rm) > 0) {
        seqCat <- seqCat[-rm]
        removeId <- ids[rm]
        ids <- ids[-rm]
        transcripts <- transcripts[-which(transcripts$transcript_id %in% 
                                              removeId)]
    }
    seqCat <- c(seqCat, stringr::str_sub(seqCat, 2), stringr::str_sub(seqCat, 
                                                                      3))
    frames <- rep(c(1, 2, 3), each = length(ids))
    ids <- c(ids, ids, ids)

    orf <- suppressWarnings(unlist(lapply(seqCat, function(x) as.character(Biostrings::translate(Biostrings::DNAString(x))))))
    orfDF <- data.frame(id = ids, aa_sequence = orf, frame = frames, 
                        stringsAsFactors = FALSE)
    orfDF$seq_length <- nchar(orfDF$aa_sequence)
    orfDF$seq_length_nt <- nchar(seqCat) + orfDF$frame - 1
    
    
    return(orfDF)
}

## 'cds_by_tx' must be a GRangesList object obtained with cdsBy(txdb, by="tx")
addCdsPhase <- function(cds_by_tx)
{
    cds_phase <- pc(rep(IntegerList(0), length(cds_by_tx)),
                    heads((3L - (cumsum(width(cds_by_tx)) %% 3L)) %% 3L, n=-1L))
    unlisted_cds_by_tx <- unlist(cds_by_tx, use.names=FALSE)
    mcols(unlisted_cds_by_tx)$cds_phase <- unlist(cds_phase, use.names=FALSE)
    relist(unlisted_cds_by_tx, cds_by_tx)
}
#produced by assigning_transcript_backbone.R
if(!exists('gtf_obj')){
    gtf_path = '/Users/annaleigh/Downloads/gencode.v38.annotation.gtf.gz'
    print(stringr::str_c(Sys.time(), " - Importing the GTF..."))
    gtf_obj =  GenomicFeatures::makeTxDbFromGFF(gtf_path,format = 'gtf')  
}
if(!exists('cds_regions')){
    cds_regions = cdsBy(gtf_obj, "tx",use.names = TRUE)
    cds_regions = unlist(cds_regions)
    
    cds_regions$transcript_id = gsub("\\..*", "", names(cds_regions))
    
}
#produced by assigning_transcript_backbone.R
with_backbone = glue::glue('{top_folder}modules_output/{experiment}.cryptic_exons.protein_coding_with_transcript_backbone.csv')
outname_peppo = glue::glue('{top_folder}modules_output/{experiment}.nmd_prediction.csv')
outname_updown = glue::glue('{top_folder}modules_output/{experiment}.cryptic_up_down.bed')

cryptic_regions_protein_coding = fread(with_backbone)
cryptic_regions_protein_coding = cryptic_regions_protein_coding %>% 
    dplyr::select(seqnames:name,transcript_id,gene_name) %>% unique() %>% 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE)



peptides = c()
with_peptide = c()
cryptic_regions_protein_coding$PTC = NA_character_
peppo_df = tibble()
all_cryptic_up_downs = GRanges()

for(i in 1:length(cryptic_regions_protein_coding)){
    if(i %% 10 == 0){
        print(i)
    }
    gene_name = cryptic_regions_protein_coding[i]$gene_name
    parent_transcript = cryptic_regions_protein_coding[i]$transcript_id
    cds_parent = cds_regions %>% filter(transcript_id == parent_transcript)
    
    if(length(cds_parent) == 0){
        next()
    }
    
    the_cryptic_exon = cryptic_regions_protein_coding[i]
    seqlevels(the_cryptic_exon) = seqlevels(cds_parent)
    
    
    new_model = GenomicRanges::reduce(sort(c(the_cryptic_exon,cds_parent)))
    
    new_model$transcript_id = parent_transcript

    new_model = GeneStructureTools::reorderExonNumbers(new_model)
    
    names(new_model) = NULL
    
    new_model = unlist(addCdsPhase(GRangesList(new_model)))
    seqlevels(new_model) = seqlevels(cds_parent)
    
    if(length( which(sort(new_model) != sort(cds_parent))) == 0){
        print(glue::glue("cryptic exon {the_cryptic_exon$cryptic_coords} in {gene_name}
                         is inside of an annotated exon in the parent transcript! Skipping"))
        next()
    }
    
    cds_parent = GeneStructureTools::reorderExonNumbers(cds_parent)
    
    cryptic_number = subsetByOverlaps(new_model,the_cryptic_exon)$exon_number
    
    cryptic_up_down = new_model %>% filter(exon_number %in% c(cryptic_number - 1, 
                                                              cryptic_number,cryptic_number + 1))
    
    cryptic_up_down = cryptic_up_down %>% mutate(type = case_when(exon_number == cryptic_number ~ 'cryptic',
                                                exon_number < cryptic_number ~ 'upstream',
                                                exon_number > cryptic_number ~ 'downstream'))
    
    cryptic_up_down$name = the_cryptic_exon$name
    
    cryptic_up_down = cryptic_up_down %>% mutate(name = glue::glue("{name}|{type}"))
    
    
    all_cryptic_up_downs = c(cryptic_up_down,all_cryptic_up_downs)
    
    normal = my_orf(cds_parent,BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38) %>% 
        dplyr::slice(1) %>% pull(aa_sequence)
    
    # tmp = GeneStructureTools::getOrfs(new_model,BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    #                                   returnLongestOnly = FALSE, 
    #                                   allFrames = TRUE,
    #                                   uORFs = TRUE)
    tmp = my_orf(new_model,BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
    tmp$gene_id = gene_name
    # tmp$nmd_prob <- notNMD::predictNMD(tmp, "prob")
    ce_aa = tmp %>% dplyr::slice(1) %>% pull(aa_sequence)
    # if(gene_name == "AARS"){
    #     break()
    # }
    if(cryptic_number == 1){
        cryptic_regions_protein_coding[i]$PTC = '5utr'
    }else if(cryptic_number == max(new_model$exon_number)){
        cryptic_regions_protein_coding[i]$PTC = '3utr'
    }else if(stringr::str_count(ce_aa,"\\*") == 1){
        both_strings = Biostrings::AAStringSet(list(Biostrings::AAString(normal), Biostrings::AAString(ce_aa)))
        aln = msa::msa(both_strings)
        ce_pos = which(strsplit(msa::msaConsensusSequence(aln)[[1]],"")[[1]] == "?")
        left_flank = (min(ce_pos) - 20) : (min(ce_pos - 1))
        left_flank = left_flank[left_flank > 0]
        right_flank = (max(ce_pos) + 1) : (max(ce_pos + 20))
        
        
        ce_peptide = paste0(strsplit(ce_aa,"")[[1]][ce_pos],collapse = "")
        
        left_flank_peptide = paste0(strsplit(ce_aa,"")[[1]][left_flank],collapse = "")
        right_flank_peptide = paste0(strsplit(ce_aa,"")[[1]][right_flank],collapse = "")
        
        if(nchar(left_flank_peptide) != 20 | nchar(right_flank_peptide) != 20){
            print(glue::glue("left and right flanks are off on {i} - {gene_name}"))
        }
        
        
        print(gene_name)
        print(ce_peptide)
        peptides = c(peptides,ce_peptide)
        with_peptide = c(with_peptide,gene_name)
        cryptic_regions_protein_coding[i]$PTC = 'no'
        the_peppos = list(gene_name,left_flank_peptide,ce_peptide,right_flank_peptide)
        names(the_peppos) = c("gene_name","left_flank_peptide","ce_peptide","right_flank_peptide")
        peppo_df = rbind(peppo_df,the_peppos)
        
    }else{
        cryptic_regions_protein_coding[i]$PTC = 'yes'
    }
    
}
beepr::beep(4)

fwrite(unique(as.data.table(cryptic_regions_protein_coding)),outname_peppo)
rtracklayer::export(unique(all_cryptic_up_downs),outname_updown)
