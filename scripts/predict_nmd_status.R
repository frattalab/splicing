library(stringi)
library(msa)
my_orf = function(transcripts, BSgenome = NULL, returnLongestOnly = TRUE, 
                  allFrames = FALSE,  
                  uORFs = FALSE) 
{
    if (allFrames == TRUE) {
        returnLongestOnly = FALSE
        longest = 1
    }
    
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
#produced by assigning_transcript_backbone.R
max_tx_protein_coding_genes
#produced by assigning_transcript_backbone.R
if(!exists(gtf_obj)){
    gtf_path = '/Users/annaleigh/Downloads/gencode.v38.annotation.gtf.gz'
    print(stringr::str_c(Sys.time(), " - Importing the GTF..."))
    gtf_obj =  GenomicFeatures::makeTxDbFromGFF(gtf_path,format = 'gtf')  
}

cds_regions = cdsBy(gtf_obj, "tx",use.names = TRUE)
cds_regions = unlist(cds_regions)

cds_regions$transcript_id = gsub("\\..*", "", names(cds_regions))


peptides = c()
with_peptide = c()
cryptic_regions_protein_coding$PTC = NA_character_
peppo_df = tibble()

for(i in 1:length(cryptic_regions_protein_coding)){
    print(i)
    gene_name = cryptic_regions_protein_coding[i]$gene_name
    print(gene_name)
    parent_transcript = cryptic_regions_protein_coding[i]$transcript_id
    cds_parent = cds_regions %>% filter(transcript_id == parent_transcript)
    
    if(length(cds_parent) == 0){
        next()
    }
    
    exon_one = cryptic_regions_protein_coding[i]
    
    new_model = sort(c(exon_one,cds_parent))
    
    new_model = GeneStructureTools::reorderExonNumbers(new_model)
    
    names(new_model) = NULL
    new_model = unlist(addCdsPhase(GRangesList(new_model)))
    
    cryptic_number = new_model %>% plyranges::filter(is.na(cds_id)) %>% 
        as.data.frame() %>% pull(exon_number)
    
    cryptic_up_down = new_model %>% filter(exon_number %in% c(cryptic_number - 1, 
                                                              cryptic_number,cryptic_number + 1))
    
    cds_parent = GeneStructureTools::reorderExonNumbers(cds_parent)
    
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
        right_flank = (max(ce_pos) + 1) : (max(ce_pos + 20))
        
        
        ce_peptide = paste0(strsplit(ce_aa,"")[[1]][ce_pos],collapse = "")
        left_flank_peptide = paste0(strsplit(ce_aa,"")[[1]][left_flank],collapse = "")
        right_flank_peptide = paste0(strsplit(ce_aa,"")[[1]][right_flank],collapse = "")
        if(nchar(left_flank_peptide) != 20 | nchar(right_flank_peptide) != 20){
            browser()
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

