#' NMDer workflow: Preparing for analysis
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE NMDer PROGRAM WORKFLOW.
#' This function will prepare output dataframes and GRanges objects for downstream workflow
#'
#' @param inputGRanges 
#' Query GRanges object
#' @param basicGRanges 
#' Reference GRanges object. 
#'
#' @return
#' Function will return a list containing:
#' (1) dataframe of a list of query transcripts for downstream workflow
#' (2) GRangesList of exon coordinates from query as dataframe
#' (3) GRangesList of exon coordinates from reference as dataframe
#' (4) GRangesList of exon coordinates of reference CDS as dataframe
#' 

prepareAnalysis <- function(inputGRanges, basicGRanges) {
  
  infoLog('Preparing databases, transcripts and output files...', logf, quiet)
  
  # prepare output dataframe by getting transcript information from inputGRanges
  report_df = inputGRanges %>% as.data.frame() %>%
    dplyr::filter(type == 'transcript') %>%
    dplyr::select(Gene_ID = gene_id, 
                  Original_Gene_ID = old_gene_id,
                  Gene_name = gene_name,
                  Match_level = match_level,
                  Transcript_ID = transcript_id,
                  seqnames,
                  start,
                  end,
                  Strand = strand) %>%
    dplyr::mutate(NMDer_ID = as.character(NA),
                  Alt_tx = as.logical(NA),
                  annotatedStart = as.logical(NA),
                  predictedStart = as.logical(NA),
                  Transcript_coord = paste0(seqnames, ':', start, '-', end))
  
  # prepare df of gencode_basic transcripts
  basicTX_df = basicGRanges %>% as.data.frame() %>%
    dplyr::filter(type == 'transcript') %>%
    dplyr::select(Gene_ID = gene_id, Ref_transcript_ID = transcript_id) %>%
    dplyr::distinct()
  
  # join report_df and basicTX_df and generate a full list of query transcripts to be compared with every reference transcript
  # and add in NMDerIDs
  combined_report_df = report_df %>%
    dplyr::left_join(basicTX_df, by = 'Gene_ID') %>%
    dplyr::select(NMDer_ID, Gene_ID, Original_Gene_ID, Match_level, 
                  Gene_Name, Transcript_ID, Ref_transcript_ID, Transcript_coord, 
                  Strand, annotatedStart, predictedStart, Alt_tx) 
  
  # prepare databases
  inputDB = makeTxDbFromGRanges(inputGRanges)
  basicDB = makeTxDbFromGRanges(basicGRanges)
  
  # prepare exon structures by transcripts/CDSs
  inputExonsbyTx = exonsBy(inputDB, by="tx", use.names=TRUE) %>% as.data.frame()
  basicExonsbyCDS = cdsBy(basicDB, by="tx", use.names=TRUE) %>% as.data.frame()
  basicExonsbyTx = exonsBy(basicDB, by="tx", use.names=TRUE) %>% as.data.frame()
  
  # remove comparisons in which the reference transcript do not have a CDS
  basicCDS = basicExonsbyCDS %>% as.data.frame() %>%
    dplyr::select(Ref_transcript_ID = group_name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(CDS = TRUE)
  
  combined_report_df = combined_report_df %>%
    dplyr::left_join(basicCDS, by = 'Ref_transcript_ID') %>%
    dplyr::filter((!is.na(CDS) & !is.na(Ref_transcript_ID)) | Match_level == 5) %>%
    dplyr::select(-CDS)
  
  # combine query transcripts with multiple comparisons to reference
  combined_report_df = combined_report_df %>% 
    group_by(Transcript_ID) %>%
    mutate(Ref_transcript_ID = paste(as.character(Ref_transcript_ID), collapse = '_')) %>%
    ungroup() %>%
    distinct(.keep_all = TRUE) %>%
    dplyr::mutate(NMDer_ID = paste0('NMDer', formatC(as.integer(row_number()), width=7, flag='0')))
  
  
  # cleanup unused variables
  rm(list = c('inputDB','basicDB'))
  return(list(combined_report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx))
}

