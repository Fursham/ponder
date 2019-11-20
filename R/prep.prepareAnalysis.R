#' _ workflow: Preparing for analysis
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW.
#' This function will prepare output dataframes and GRanges objects for the _Program 
#'
#' @param inputGRanges 
#' GRanges object containing transcript and exon features. It is recommended to create this object
#' by importing gtf/gff3 files using rtracklayer::import
#' @param basicGRanges 
#' #' GRanges object containing reference transcript,exon and CDS features. 
#' It is recommended to create this object by importing gtf/gff3 files using rtracklayer::import
#'
#' @return
#' Function will return a list containing:
#' (1) dataframe of a list of assembled transcripts and its information for output
#' (2) dataframe of a list of reference CDS transcripts
#' (3) GRangesList of exon coordinates of assembled transcripts grouped by transcript name
#' (4) GRangesList of exon coordinates of reference CDS grouped by transcript name
#' 
#' @export
#' @import dplyr
#' @import GenomicFeatures
#'
#' @examples
prepareAnalysis <- function(inputGRanges, basicGRanges) {
  
  infoLog('Preparing databases, transcripts and output files...', logf, quiet)
  
  # prepare output dataframe by getting transcript information from inputGRanges
  report_df = inputGRanges %>% as.data.frame() %>%
    dplyr::filter(type == 'transcript') %>%
    dplyr::select(Gene_ID = gene_id, 
                  Original_Gene_ID = old_gene_id,
                  Gene_Name = gene_name,
                  Match_level = match_level,
                  Transcript_ID = transcript_id,
                  seqnames,
                  start,
                  end,
                  Strand = strand) %>%
    dplyr::mutate(NMDer_ID = as.character(NA),
                  Alt_tx = as.logical(NA),
                  annotatedStart = as.logical(NA),
                  predictedStart = as.logical(NA)) %>%
    tidyr::unite(Tx_coord., c(start, end), sep = '-') %>%
    tidyr::unite(Tx_coord, c(seqnames, Tx_coord.), sep = ':') 
  
  # prepare df of gencode_basic transcripts
  basicTX_df = basicGRanges %>% as.data.frame() %>%
    dplyr::filter(type == 'transcript') %>%
    dplyr::select(Gene_ID = gene_id, Ref_TX_ID = transcript_id) %>%
    dplyr::distinct()
  
  # join report_df and basicTX_df and generate a full list of query transcripts to be compared with every reference transcript
  # and add in NMDerIDs
  combined_report_df = suppressMessages(report_df %>%
                                          dplyr::left_join(basicTX_df) %>%
                                          dplyr::select(NMDer_ID, Gene_ID, Original_Gene_ID, Match_level, 
                                                        Gene_Name, Transcript_ID, Ref_TX_ID, Tx_coord, 
                                                        Strand, annotatedStart, predictedStart, Alt_tx)) 
  
  # prepare databases
  inputDB = makeTxDbFromGRanges(inputGRanges)
  basicDB = makeTxDbFromGRanges(basicGRanges)
  
  # prepare exon structures by transcripts/CDSs
  inputExonsbyTx = exonsBy(inputDB, by="tx", use.names=TRUE) %>% as.data.frame()
  basicExonsbyCDS = cdsBy(basicDB, by="tx", use.names=TRUE) %>% as.data.frame()
  basicExonsbyTx = exonsBy(basicDB, by="tx", use.names=TRUE) %>% as.data.frame()
  
  # remove comparisons in which the reference transcript do not have a CDS
  basicCDS = basicExonsbyCDS %>% as.data.frame() %>%
    dplyr::select(Ref_TX_ID = group_name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(CDS = TRUE)
  
  combined_report_df = suppressMessages(combined_report_df %>%
                                          dplyr::left_join(basicCDS) %>%
                                          dplyr::filter((!is.na(CDS) & !is.na(Ref_TX_ID)) | Match_level == 5) %>%
                                          dplyr::select(-CDS))
  
  # combine query transcripts with multiple comparisons to reference
  combined_report_df = combined_report_df %>% 
    group_by(Transcript_ID) %>%
    mutate(Ref_TX_ID = paste(as.character(Ref_TX_ID), collapse = '_')) %>%
    ungroup() %>%
    distinct(.keep_all = TRUE) %>%
    dplyr::mutate(NMDer_ID = paste0('NMDer', formatC(as.integer(row_number()), width=7, flag='0')))
  
  
  # cleanup unused variables
  rm(list = c('inputDB','basicDB'))
  return(list(combined_report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx))
}

