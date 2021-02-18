#!/usr/bin/env Rscript
library(GenomicFeatures)
library(tidyverse)
library(data.table)
library("optparse")
# use a function from splicejam to find the exons which end on any of the junctions
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
    splicejamobject = splicejam::closestExonToJunctions(
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

combine_exons_junctions = function(transcripts, junctions, output_name){
  ####read in the scallop GTF
  gtf = transcripts
  # gtf = "/Users/annaleigh/Downloads/gencode.v37.basic.annotation.gtf"
  scallop_form =  GenomicFeatures::makeTxDbFromGFF(gtf,format = 'gtf')
  ####read in the bed file
  bed_junctions = data.table::fread(junctions,skip = 1)
  bed_junctions = bed_junctions %>%
      mutate(control_psi = as.numeric(stringr::str_split(V4,"_0",n = 2, simplify = T)[,2])) %>%
      mutate(control_psi = ifelse(is.na(control_psi),0,control_psi)) %>%
      separate(V4,into = c('gene_name','type'),sep = ":") %>%
      separate(type, into = "type",sep = "_0")
  #### get the junctions
  # tx_jnc = as.data.table(intronsByTranscript(scallop_form,use.names=T))
  ### get the exons
  exon_jnc = as.data.table(exonsBy(scallop_form,use.names=T,by = 'tx'))
  ### add an exon column
  exon_jnc[,paste_into_igv_exon := paste0(seqnames, ":",start, "-",end)]
  exon_jnc[,exon_name := paste0(group_name,":",exon_rank)]

  ##flatten the exons
  scallop_exons_gr = makeGRangesFromDataFrame(exon_jnc,keep.extra.columns = T)
  bed_gf = makeGRangesFromDataFrame(bed_junctions,
                                    start.field = 'V2',
                                    end.field = 'V3',
                                    seqnames.field = 'V1',
                                    strand.field = 'V6',
                                    keep.extra.columns = T)

  exons_on_junction_ends = build_all_potential_exon_table(bed_gf,scallop_exons_gr)

  exons_on_junction_ends = annotate_all_possible_exon_table(exons_on_junction_ends,exon_jnc)

  exons_on_junction_ends = exons_on_junction_ends %>%
      left_join(bed_junctions, by = c("seqnames" = "V1",
                                      "start" = "V2",
                                      "end" = "V3")) %>%
      unique()

  exons_on_junction_ends = exons_on_junction_ends %>% mutate(cryptic_exon = case_when(strand == "+" & type == "novel_acceptor" ~ Exon_target,
                                           strand == "+" & type == "novel_donor" ~ Exon_source,
                                           strand == "-" & type == "novel_acceptor" ~ Exon_source,
                                           strand == "-" & type == "novel_donor" ~ Exon_target))



  cryptic_exons = unique(exons_on_junction_ends[grepl("novel",type) & V5 > 0.1 & control_psi < 0.05][,.(cryptic_exon,gene_name)])



  cryptic_exons = unique(exons_on_junction_ends[grepl("novel",type) & V5 > 0.1 &
                                    control_psi < 0.05][,.(cryptic_exon,gene_name,V5)][,tdpKD_psi := mean(V5),by = .(cryptic_exon,gene_name)])


  cryptic_exons %>%
      mutate(strand = str_sub(cryptic_exon,-1,-1)) %>%
      separate(cryptic_exon, into = c("seqname","start","end")) %>%
      filter(!is.na(start)) %>%
      dplyr::select(-V5) %>%
      unique() %>%
      fwrite(output_name,col.names = F, sep = "\t")
}
####end of the helper functions####
option_list = list(
    make_option(c("-t", "--transcripts"), type="character", default=NULL,
                help="assembled transcripts from the transcriptome_assembly workflow", metavar="character"),
    make_option(c("-d", "--deltabed"), type="character", default=NULL,
                help="one of the deltaPSI beds produced by this workflow", metavar="character"),
    make_option(c("-o", "--outputname"), type="character",
                help="output file name - no extension 2 files will be written", metavar="character",default = "cryptic_exons.bed")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

combine_exons_junctions(transcripts = opt$transcripts, junctions = opt$deltabed, output_name = opt$out)
