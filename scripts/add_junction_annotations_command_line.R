#!/usr/bin/env Rscript

library("optparse")
library(GenomicRanges)
####Start of the helper functions####
#' Annotate junctions using reference annotation
#'
#' \code{annotate_junc_ref} annotates junctions by whether their start and end
#' position precisely overlaps with a known exon boundary. Using this
#' information along with the strand, junctions are categorised into
#' "annotated", "novel_acceptor", "novel_donor", "novel_combo",
#' "novel_exon_skip", "ambig_gene" and "none".
#'
#' @param junc_metadata junction metadata in a
#'   \code{\link[GenomicRanges]{GRanges}} format, the essential component being
#'   the junction co-ordinates.
#' @param gtf either path to gtf or object of class \code{ensemblGenome} loaded
#'   using \code{refGenome}.
#'
#' @return junction metadata as a \code{\link[GenomicRanges]{GRanges}} object
#'   with additional columns that detail overlapping genes/transcripts/exons and
#'   junction categories.
#'
#' @export
annotate_junc_ref <- function(junc_metadata, gtf){
    
    ##### Check user input is correct #####
    
    if(!("GRanges" %in% class(junc_metadata))) stop("junction_metadata must be in a GRanges format")
    
    if(all(!(c("character", "TxDb") %in% class(gtf)))){
        
        stop("gtf must either be a path to the .gtf file or a pre-loaded gtf of class TxDb")
        
    }
    
    ##### Extract annotated exons/junctions co-ordinates from gtf #####
    
    print(stringr::str_c(Sys.time(), " - Obtaining co-ordinates of annotated exons and junctions from gtf..."))
    
    if(class(gtf) == "character"){
        
        print(stringr::str_c(Sys.time(), " - Importing gtf..."))
        
        # import gtf using refGenome, needed to obtain the annotated splice junctions easily
        # ref <- refGenome::ensemblGenome()
        # refGenome::basedir(ref) <- dirname(gtf)
        # refGenome::read.gtf(ref, gtf %>% stringr::str_replace(".*/", ""))
        ref =  GenomicFeatures::makeTxDbFromGFF(gtf,format = 'gtf')
        
    }else if(class(gtf) == "TxDb"){
        
        ref <- gtf
        
    }
    
    # ref_exons <- ref@ev$gtf[ref@ev$gtf$feature == "exon",] %>% GRanges()
    # ref_junc <- refGenome::getSpliceTable(ref)
    # ref_junc <- ref_junc@ev$gtf
    ref_exons <- exonsBy(gtf,'tx',use.name = TRUE) ##get the exons by the transcript
    
    ref_junc <- BiocGenerics::setdiff(range(ref_exons), ref_exons) ### junctions are just teh range of eeach transcript minus where the exons are
    ##### Obtain annotation through overlapping exons #####
    print(stringr::str_c(Sys.time(), " - Adding gene and transcript information to exons..."))
    ref_exons <- .add_gene_information(ref_exons, gtf)
    
    print(stringr::str_c(Sys.time(), " - Getting junction annotation using overlapping exons..."))
    
    junc_metadata <- .get_ref_exons_annot(junc_metadata, ref_exons)
    
    ##### Tidy annotation - collapse gene annotation to per junction and infer strand #####
    
    print(stringr::str_c(Sys.time(), " - Tidying junction annotation..."))
    
    junc_metadata <- .tidy_junc_annot(junc_metadata)
    
    ##### Derive junction categories using strand & overlapping exon annotation #####
    
    print(stringr::str_c(Sys.time(), " - Deriving junction categories..."))
    
    junc_metadata <- .classify_junc(junc_metadata, ref_junc)
    
    print(stringr::str_c(Sys.time(), " - done!"))
    
    return(junc_metadata)
    
}

#' Extracts annotation from the reference gtf
#'
#' \code{.convert_hits_to_annot_ref} will find the overlap of each junction with
#' the annotated exons. Then, for each corresponding hit, annotates each
#' junction with the strand/exon/transcript/gene from the reference annotation.
#'
#' @inheritParams annotate_junc_ref
#'
#' @param ref_exons annotated exons imported from the gtf.
#' @param junc_start_end start and end of the junctions returned from
#'   \code{\link{.get_gr_for_start_end}}.
#' @param ref_exons_start_end start and end of the exons returned from
#'   \code{\link{.get_gr_for_start_end}}.
#' @param ref_cols column names matching the annotation columns from the
#'   reference.
#'
#' @return junction metadata with annotation.
.get_ref_exons_annot <- function(junc_metadata,
                                 ref_exons,
                                 ref_cols = c("strand", "transcript_id", "exon_name","exon_rank","gene_id")){
    
    # match junctions to exon definitions
    start(junc_metadata) <- start(junc_metadata) - 1
    end(junc_metadata) <- end(junc_metadata) + 1
    
    # make a gr where each junc/exon is marked by only a start or end co-ordinate
    junc_start_end <- .get_gr_for_start_end(junc_metadata)
    ref_exons_start_end <- .get_gr_for_start_end(ref_exons)
    
    for(start_end in c("start", "end")){
        
        # only get hits between junc start/exon end or junc end/exon start
        # the other way (e.g. junc end/exon end) should not happen (only 0.05% of the data)
        end_start <- ifelse(start_end == "start", "end", "start")
        
        # avoid seqlevel non-overlap warnings
        suppressWarnings(
            junc_exon_hits <- findOverlaps(query = junc_start_end[[start_end]],
                                           subject = ref_exons_start_end[[end_start]],
                                           type = "equal",
                                           ignore.strand = F)
        )
        
        # set junc_hits to factor with levels containing all junction indexes
        # so split(drop = F) keeps all junctions
        # not only those which precisely overlap an exon boundary
        junc_hits_fct <- queryHits(junc_exon_hits) %>%
            factor(levels = 1:length(junc_metadata))
        
        for(j in seq_along(ref_cols)){
            
            # extract the values from exon metadata column of interest
            if(ref_cols[j] != "strand"){
                
                ref_col_values <- ref_exons %>%
                    mcols() %>%
                    .[[ref_cols[j]]]
                
            }else{
                
                # if strand extract strand
                ref_col_values <- strand(ref_exons)
                
            }
            
            mcols(junc_metadata)[stringr::str_c(ref_cols[j], "_", start_end)] <- ref_col_values %>%
                .[subjectHits(junc_exon_hits)] %>% # subset the exons by those that overlap juncs
                split(junc_hits_fct, drop = F) %>% # split into groups based on index of overlapping junc
                CharacterList() %>%
                unique() # parrallelised unique - remove duplicates when for example strand if junc overlaps >1 exon
            
        }
        
    }
    
    # convert junc co-ords back to intron definitions
    start(junc_metadata) <- start(junc_metadata) + 1
    end(junc_metadata) <- end(junc_metadata) - 1
    
    return(junc_metadata)
    
}

#' Add the gene name/gene_id onto the transcript annotated object
#'
#' \code{.add_gene_information} will add the gene name and id onto the annotated
#' transcript/exon the overlap of each junction with
#'
#' @inheritParams annotate_junc_ref
#'
#' @param gtf reference gtf

#' @return junction metadata with gene_id annotation
.add_gene_information <- function(ref_exons,
                                  gtf){
    
    #unlist the exons
    ref_exons = unlist(ref_exons)
    ref_exons$transcript_id = names(ref_exons) #add on the transcript_id to the flattened list
    #unlist the exons
    tx_by_gene = unlist(transcriptsBy(gtf,'gene'))
    tx_by_gene$gene_id = names(tx_by_gene)
    tx_by_gene$transcript_id  = tx_by_gene$tx_name
    tx_by_gene$tx_name  = NULL
    tx_by_gene$tx_id  = NULL
    
    
    mcols(ref_exons) = mcols(ref_exons) %>% as.data.frame()%>%
        dplyr::left_join(as.data.frame(mcols(tx_by_gene)),by = 'transcript_id') %>% unique()
    
    
    return(ref_exons)
    
}



#' Tidy junction annotation
#'
#' \code{.tidy_junc_annot} merges the gene and strand details from the start and
#' end into one column per junction. Then, combines strand information from the
#' original rna-seq based and that from overlapping annotation.
#'
#' @inheritParams annotate_junc_ref
#'
#' @param cols_to_merge
#'
#' @return
.tidy_junc_annot <- function(junc_metadata, cols_to_merge = c("strand", "gene_id")){
    
    # collapse gene/strand columns to per junc instead of per start/end for easier querying
    for(col in cols_to_merge){
        
        mcols(junc_metadata)[[stringr::str_c(col, "_junc")]] <-
            .merge_lists(mcols(junc_metadata)[[stringr::str_c(col, "_start")]],
                         mcols(junc_metadata)[[stringr::str_c(col, "_end")]])
        
    }
    
    # replacing empty strands ("none") and those with >1 strand ("ambig_gene") with "*"
    # this ensures each vector in CharacterList is of length 1
    # so can be unlisted and length(chr_list) == length(unlist(chr_list))
    mcols(junc_metadata)[["strand_junc"]][lengths(mcols(junc_metadata)[["strand_junc"]]) == 0] <- "*"
    mcols(junc_metadata)[["strand_junc"]][lengths(mcols(junc_metadata)[["strand_junc"]]) > 1] <- "*"
    strand_annot <- unlist(mcols(junc_metadata)[["strand_junc"]])
    
    # compare the strand obtained from annotation strand to original strand
    orig_strand <- as.character(strand(junc_metadata))
    
    # salvage situations when either original or annotation strand is "*" and the other is "+" or "-"
    strand(junc_metadata) <- dplyr::case_when(orig_strand == strand_annot ~ orig_strand,
                                              orig_strand == "*" & strand_annot != "*" ~ strand_annot,
                                              strand_annot == "*" & orig_strand != "*" ~ orig_strand)
    
    # remove to avoid confusion between strand() and strand_junc
    mcols(junc_metadata)[["strand_junc"]] <- NULL
    
    return(junc_metadata)
    
}

#' Classifys junctions
#'
#' \code{.classify_junc} categories junctions into "annotated",
#' "novel_acceptor", "novel_donor", "novel_combo", "exon_skip", "ambig_gene" and
#' "none" using information from annotation and strand. Adds two additional
#' columns
#'
#' @inheritParams annotate_junc_ref
#'
#' @return junction metadata with additional columns informing whether that
#'   junction is present within annotation and it's category.
.classify_junc <- function(junc_metadata, ref_junc){
    
    ref_junc = unlist(ref_junc)
    ref_junc$transcript_id = names(ref_junc)
    # find whether junction is found in splice table
    
    # ref_junc_gr <- ref_junc %>%
    #   as.data.table()
    #   dplyr::mutate(start = start + 1, # match exon boundaries to intron co-ords
    #                 end = end -1)
    #   GRanges() %>%
    #   unique()
    
    # avoid diff seqlevels warning
    suppressWarnings(annot_hits <- findOverlaps(query = junc_metadata,
                                                subject = ref_junc,
                                                type = "equal"))
    
    mcols(junc_metadata)[["junc_in_ref"]] <- 1:length(junc_metadata) %in% queryHits(annot_hits)
    
    # classify junctions
    # separate strand out for readability
    strand_junc <- as.character(strand(junc_metadata))
    mcols(junc_metadata)[["junc_cat"]] <-
        dplyr::case_when(mcols(junc_metadata)[["junc_in_ref"]] == T ~ "annotated",
                         lengths(mcols(junc_metadata)[["gene_id_junc"]]) == 0 ~ "none",
                         lengths(mcols(junc_metadata)[["gene_id_junc"]]) > 1 ~ "ambig_gene", # after these checks lengths(gene_id_junc) must equal 1
                         lengths(mcols(junc_metadata)[["gene_id_start"]]) > 0 & lengths(mcols(junc_metadata)[["gene_id_end"]]) > 0 ~ "novel_combo",
                         strand_junc == "+" & lengths(mcols(junc_metadata)[["gene_id_start"]]) > 0 ~ "novel_acceptor",
                         strand_junc == "-" & lengths(mcols(junc_metadata)[["gene_id_start"]]) > 0 ~ "novel_donor",
                         strand_junc == "+" & lengths(mcols(junc_metadata)[["gene_id_end"]]) > 0 ~ "novel_donor",
                         strand_junc == "-" & lengths(mcols(junc_metadata)[["gene_id_end"]]) > 0 ~ "novel_acceptor")
    
    # split the novel_combo into novel_combo and novel_exon_skip
    # do this separately to save time - case_when evaluates each condition across all junctions
    # since each column is called as mcols(junc_metadata)[["col"]]
    # converting this data_frame() parses CharacterLists to list(), so also non-optimal solution
    mcols(junc_metadata)[["index_tmp"]] <- 1:length(junc_metadata)
    novel_combo <- junc_metadata[mcols(junc_metadata)[["junc_cat"]] == "novel_combo"]
    exon_skip_indexes <- mcols(novel_combo)[["index_tmp"]][any(novel_combo$transcript_id_start %in% novel_combo$transcript_id_end)]
    mcols(junc_metadata)[["junc_cat"]][exon_skip_indexes] <- "novel_exon_skip"
    mcols(junc_metadata)[["index_tmp"]] <- NULL
    
    return(junc_metadata)
    
}



####end of the helper functions####
option_list = list(
    make_option(c("-d", "--deltapsi"), type="character", default=NULL, 
                help="parsed deltaPSI table produced by parse_voila_delta_tsv script", metavar="character"),
    make_option(c("-o", "--out"), type="character", default=".", 
                help="output file name", metavar="character"),
    make_option(c("-g", "--gtf"), type="character", default=NULL, 
                help="the GTF to annotate to", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

annotate_junctions <- function(parsed_file,output_filepath,gtf){
    ####tessst####
    parsed_splicing = data.table::fread(parsed_file)
    txdb_gtf =  GenomicFeatures::makeTxDbFromGFF(gtf,format = 'gtf') #do we need this still???
    
    # use this superset of all the exons either identified by Whippets
    # or identified by Scallop
    # use a function from splicejam to find the exons which end on any of the junctions
    parsed_granges = makeGRangesFromDataFrame(ward_splicing,
                                            start.field = "junc_start",
                                            end.field = "junc_end",
                                            keep.extra.columns = TRUE)
    
    start(parsed_granges) <- start(parsed_granges) + 1 #the annotate junc ref tool is based on STAR SJ's so there's an internal -1, +1 that happens
    end(parsed_granges) <- end(parsed_granges) - 1
    
    annotated_junctions = annotate_junc_ref(parsed_granges,txdb_gtf)
    #####put the coordinates back
    start(annotated_junctions) <- start(annotated_junctions) - 1 #the annotate junc ref tool is based on STAR SJ's so there's an internal -1, +1 that happens
    end(annotated_junctions) <- end(annotated_junctions) + 1
    


}