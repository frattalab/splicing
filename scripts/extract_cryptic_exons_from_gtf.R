#!/usr/bin/env Rscript
library(GenomicFeatures)
library(tidyverse)
library(data.table)
library("optparse")
# use a function from splicejam to find the exons which end on any of the junctions

source("scripts/splicejam_closestExonToJunctions.R")

build_all_potential_exon_table = function(grange_object,exon_annotation){

    ##check that exon_coords exists as a metadata column




    if(!("exon_coords" %in% (colnames(mcols(exon_annotation))))){
        exon_annotation$exon_coords = paste0(as.character(seqnames(exon_annotation)),
                                             ":",
                                             as.character(start(exon_annotation)),
                                             "-",
                                             as.character(end(exon_annotation)),
                                             ":",
                                             as.character(strand(exon_annotation)))
    }
    ##check that the exons have names, splice jam needs it
    if(is.null(names(exon_annotation))){
        names(exon_annotation) = exon_annotation$exon_coords
    }
    ####build a table of all potential exons for all junctions####
    splicejamobject = splicejam_closestExonToJunctions(
        grange_object,
        exon_annotation,
        flipNegativeStrand = TRUE,
        sampleColname = "sample_id",
    )

    ####build a table of all potential exons for all junctions####
    spliceStartExonEndD = splicejamobject$spliceStartExonEndD
    spliceEndExonStartD = splicejamobject$spliceEndExonStartD

    #retrieve the junctions that we queried, the splice start
    inputSpliceStarts = as.data.table(grange_object[spliceStartExonEndD$queryHits,]) %>%
        dplyr::select(seqnames,
                      start,
                      end,
                      strand)

    # get all the exons taht end on those splice starts
    inputExonsEnds = as.data.table(exon_annotation[spliceStartExonEndD$subjectHits,]) %>%
        dplyr::select(exon_coords)

    allPossibleSourceExons = cbind(inputSpliceStarts,inputExonsEnds,spliceStartExonEndD$strandedDistance) %>%
        dplyr::rename(Exon_source = exon_coords) %>%
        filter(V3 == 0) %>%
        dplyr::select(-V3) %>%

        unique()

    ####now for the TargetExons
    ####all possible target exons ####
    inputSpliceEnds = as.data.table(grange_object[spliceEndExonStartD$queryHits,]) %>%
        dplyr::select(seqnames,
                      start,
                      end,
                      strand)

    inputExonsStarts = as.data.table(exon_annotation[spliceEndExonStartD$subjectHits,]) %>%
        dplyr::select(exon_coords)

    allPossibleTargetExons = cbind(inputSpliceEnds,inputExonsStarts,spliceEndExonStartD$strandedDistance) %>%
        dplyr::rename(Exon_target = exon_coords) %>%

        filter(V3 == 0) %>% #only take distance zero
        dplyr::select(-V3) %>%
        unique()

    ###merge together the source and target
    allPossibleJunctionExon = allPossibleSourceExons %>%
        left_join(allPossibleTargetExons) %>%
        filter(!is.na(Exon_target) & !is.na(Exon_source)) %>%
        unique() %>% as.data.table()

    allPossibleJunctionExon[,paste_into_igv_junction := paste0(seqnames, ":",start, "-",end)]

    return(allPossibleJunctionExon)
}


annotate_all_possible_exon_table = function(possible_exons, exon_by_transcript_dt){
    if(!("exon_coords" %in% colnames(exon_by_transcript_dt))){
        exon_by_transcript_dt[,exon_coords := paste0(seqnames,":",start,"-",end,":",strand)]
    }


    annotated_pe = possible_exons %>%
        left_join(exon_by_transcript_dt %>% dplyr::select(exon_coords,
                                                          group_name,
                                                          exon_rank) %>% unique(),
                  by = c("Exon_source" = "exon_coords")) %>%
        dplyr::rename(transcript_source = group_name,
                      exon_rank_source = exon_rank) %>%
        left_join(exon_by_transcript_dt %>% dplyr::select(exon_coords,
                                                          exon_rank) %>% unique(),
                  by = c("Exon_target" = "exon_coords")) %>%
        dplyr::rename(exon_rank_target = exon_rank) %>% unique()

    annotated_pe = annotated_pe %>% mutate( paste_into_igv_junction := paste0(seqnames, ":",start, "-",end))
    return(annotated_pe)
}

combine_exons_junctions = function(transcripts, junctions, output_name, contrast_min = 0.15, base_max = 0.05){
  ### read in the parsed annotated junctions
  delta_psi = fread(junctions)
  delta_psi = delta_psi %>%
      select(seqnames:lsv_type,junc_cat) %>% 
      select(-lsv_type)
  
  cryptic_junctions = delta_psi %>% 
      filter(delta_psi[,12] < base_max & delta_psi[,13] > contrast_min)
  ####read in the transcriptome assembled GTF
  # gtf = "/Users/annaleigh/Downloads/gencode.v37.basic.annotation.gtf"
  gtf_txdb =  GenomicFeatures::makeTxDbFromGFF(transcripts,format = 'gtf')
  #### get the junctions
  # tx_jnc = as.data.table(intronsByTranscript(gtf_txdb,use.names=T))
  ### get the exons
  exon_jnc = as.data.table(exonsBy(gtf_txdb,use.names=T,by = 'tx'))
  ### add an exon column
  exon_jnc[,paste_into_igv_exon := paste0(seqnames, ":",start, "-",end)]
  exon_jnc[,exon_name := paste0(group_name,":",exon_rank)]
  exon_jnc = exon_jnc %>% 
      group_by(group_name) %>% 
      mutate(last_exon = (exon_rank == max(exon_rank))) %>% 
      mutate(first_exon = (exon_rank == min(exon_rank))) %>% 
      as.data.table()
  ##flatten the exons
  flattened_exons = makeGRangesFromDataFrame(exon_jnc,keep.extra.columns = T)
  ##turn the junctions into a grange the exons
  
  cj_gr = makeGRangesFromDataFrame(cryptic_junctions,
                                    keep.extra.columns = T)

  exons_on_junction_ends = build_all_potential_exon_table(cj_gr,flattened_exons)

  exons_on_junction_ends = annotate_all_possible_exon_table(exons_on_junction_ends,exon_jnc)

  exons_on_junction_ends = exons_on_junction_ends %>%
      left_join(cryptic_junctions, by = c("seqnames",
                                      "start",
                                      "end", "strand")) %>%
      unique()

  exons_on_junction_ends = exons_on_junction_ends %>% 
      mutate(cryptic_exon = case_when(strand == "+" & junc_cat == "novel_acceptor" ~ Exon_target,
                                           strand == "+" & junc_cat == "novel_donor" ~ Exon_source,
                                           strand == "-" & junc_cat == "novel_acceptor" ~ Exon_source,
                                           strand == "-" & junc_cat == "novel_donor" ~ Exon_target))


  cryptic_exons = exons_on_junction_ends %>% 
      filter(!is.na(cryptic_exon)) %>% 
      left_join(exon_jnc[,.(last_exon,first_exon,exon_coords)],by = c("cryptic_exon" = "exon_coords")) %>% 
      select(gene_name,paste_into_igv_junction,probability_changing:first_exon) %>% 
      unique()
  
  
  colnames(cryptic_exons) = c("gene_name", "paste_into_igv_junction", "probability_changing", 
    "probability_non_changing", "baseline_PSI", "contrast_PSI", 
    "junc_cat", "cryptic_exon", "last_exon", "first_exon")
  
  
  cryptic_exons = cryptic_exons %>% 
      group_by(cryptic_exon) %>% 
      mutate(baseline_avg_PSI = mean(baseline_PSI),
             contrast_avg_PSI = mean(contrast_PSI)) %>% 
      ungroup() %>% 
      as.data.table()


  cryptic_exons %>%
      mutate(strand = str_sub(cryptic_exon,-1,-1)) %>%
      separate(cryptic_exon, into = c("seqname","start","end")) %>%
      filter(!is.na(start)) %>%
      mutate(name = glue::glue("{baseline_avg_PSI}|{contrast_avg_PSI}|{gene_name}")) %>% 
      rtracklayer::export(glue::glue("{outputname}.bed"))
  
  cryptic_exons %>% 
      fwrite(glue::glue("{outputname}.csv"))
}

####end of the helper functions####
option_list = list(
    make_option(c("-t", "--transcripts"), type="character", default=NULL,
                help="assembled transcripts from the transcriptome_assembly workflow", metavar="character"),
    make_option(c("-d", "--delta"), type="character", default=NULL,
                help="one of the deltaPSI annotated CSV files produced by this workflow", metavar="character"),
    make_option(c("-o", "--outputname"), type="character",
                help="output file name - no extension - 2 files will be written", metavar="character",default = "cryptic_exons.bed")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

combine_exons_junctions(transcripts = opt$transcripts, junctions = opt$delta, output_name = opt$outputname)

assign("temp_bed",bed_gf)
start(temp_bed) = start(temp_bed) + 1
# end(temp_bed) = end(temp_bed) + 1
exons_on_junction_ends = build_all_potential_exon_table(temp_bed,flattened_exons)
nrow(exons_on_junction_ends)
