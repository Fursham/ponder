#' Title
#'
#' @param queryID 
#' @param refID 
#' @param gene_id 
#' @param NMDer_ID 
#' @param inputExonsbyTx 
#' @param basicExonsbyTx 
#' @param basicExonsbyCDS 
#'
#' @return
#'
#' @examples
getBestRef <- function(queryID, refID, gene_id, NMDer_ID, inputExonsbyTx, basicExonsbyTx, basicExonsbyCDS) {
  
  # prepare output list
  output = list('queryGRanges' = NA,
                'basicTxGRanges' = NA,
                'basicCDSGRanges' = NA,
                'report' = list(
                  'Ref_transcript_ID' = as.character(NA),
                  'Coverage' = 0,
                  'ORF_considered' = NA,
                  'ORF_start' = 'Not found',
                  'ORF_found' = FALSE))
  
  refID = strsplit(refID, split = '_')  #split string of ID into a list
  
  # prepare GRanges for intersection
  queryGRanges = inputExonsbyTx %>% as.data.frame() %>%
    dplyr::filter(group_name == queryID) %>%
    dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start)) %>% 
    GenomicRanges::makeGRangesListFromDataFrame(keep.extra.columns = TRUE, split = 'group_name')
  
  basicTxGRanges = basicExonsbyTx %>% 
    dplyr::filter(group_name %in% refID[[1]]) %>%
    dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start))  %>% 
    GenomicRanges::makeGRangesListFromDataFrame(keep.extra.columns = TRUE, split = 'group_name')
  
  # intersect query to the list of reference
  overlapHits = IRanges::mergeByOverlaps(queryGRanges, basicTxGRanges)
  
  # mapply function to return Coverage value
  overlapHitsMeta = base::mapply(function(x,y){
    widthQuery = sum(IRanges::width(x))
    widthRef = sum(IRanges::width(y))
    aveWidth = (widthQuery + widthRef) / 2
    commonCoverage = sum(IRanges::width(GenomicRanges::reduce(GenomicRanges::intersect(x, y))))
    Shared_coverage = commonCoverage/aveWidth
    return(Shared_coverage)
  }, overlapHits$queryGRanges, overlapHits$basicTxGRanges) %>%
    as.data.frame() # output list as a dataframe
  
  # return
  if(nrow(overlapHitsMeta) == 0) {
    return(output)
  }
  
  # append reference tx IDs, sort dataframe and select best reference
  overlapHitsMeta$Ref_transcript_ID = names(overlapHits$basicTxGRanges)
  names(overlapHitsMeta) = c('Coverage', 'Ref_transcript_ID')
  overlapHitsMeta = overlapHitsMeta %>%
    dplyr::arrange(desc(Coverage)) %>%
    dplyr::select(Ref_transcript_ID, Coverage)
  
  outBestRef = overlapHitsMeta[1,]
  
  # prepare GRanges for output
  output$queryGRanges = inputExonsbyTx %>% as.data.frame() %>%
    dplyr::filter(group_name == queryID) %>%
    dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start)) %>% 
    dplyr::mutate(transcript_id = NMDer_ID, gene_id = gene_id) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  output$basicCDSGRanges = basicExonsbyCDS %>% 
    dplyr::filter(group_name %in% outBestRef$Ref_transcript_ID) %>%
    dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start)) %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  output$basicTxGRanges = basicExonsbyTx %>% 
    dplyr::filter(group_name %in% outBestRef$Ref_transcript_ID) %>%
    dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start))  %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  # update reference transcript and coverage
  output$report$Ref_transcript_ID = outBestRef$Ref_transcript_ID
  output$report$Coverage = outBestRef$Coverage
  
  # set query ORF if it is similar to reference
  if(outBestRef$Coverage == 1) {
    # reformat CDSGRanges
    newbasicCDSGRanges = output$basicCDSGRanges %>% as.data.frame() %>%
      dplyr::mutate(type = 'CDS', 
                    gene_id = gene_id, 
                    transcript_id = NMDer_ID) %>%
      dplyr::mutate(phase = cumsum(width%%3)%%3) %>%
      dplyr::select(seqnames:end, strand, type, phase, gene_id, transcript_id)
    newbasicCDSGRanges$phase = c(0, head(newbasicCDSGRanges$phase, - 1))
    newbasicCDSGRanges = GenomicRanges::makeGRangesFromDataFrame(newbasicCDSGRanges, keep.extra.columns = TRUE)
    
    # update output variables
    output$report$ORF_considered = newbasicCDSGRanges
    output$report$ORF_start = 'Annotated'
    output$report$ORF_found = TRUE
  }
  return(output)
}