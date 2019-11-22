#' _ workflow: Test transcript for NMD feature against a reference CDS
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW.
#' This function will compare the query transcript with a reference CDS and tests
#' if insertion of alternative segments into the CDS generates an NMD substrate
#'
#' @param knownCDS 
#' @param queryTx 
#' @param refsequence 
#' @param PTC_dist 
#' @param nonClassicalNMD 
#'
#' @return
#'
#' @examples
getORF <- function(knownCDS, queryTx, refsequence, gene_id, transcript_id) {
  
  # prep output list
  output = list(ORF_considered = as.character(NA),
                Alt_tx  = as.logical(NA),
                annotatedStart = as.logical(NA),
                predictedStart = as.logical(NA))
  
  # precheck for annotated start codon on query transcript and update output
  pre_report = testTXforStart(queryTx, knownCDS, full.output=TRUE)
  output = utils::modifyList(output, pre_report["annotatedStart"])
  
  # return if there is no shared exons between transcript and CDS
  if (is.na(pre_report$txrevise_out[1])) {
    return(output)
  } 
  
  # attempt to reconstruct CDS for transcripts with unannotated start
  if ((pre_report$annotatedStart == FALSE) |
      (pre_report$annotatedStart == TRUE & pre_report$firstexlen < 3)) {
    
    diffSegments = pre_report$txrevise_out
    queryStrand = as.character(strand(knownCDS))[1]
    pre_report$ORF = sort(append(
      diffSegments$shared_exons, c(
        reduce(diffSegments[[1]][diffSegments[[1]]$upstream != TRUE]),
        reduce(diffSegments[[2]][diffSegments[[2]]$upstream == TRUE]))),
      decreasing = queryStrand == '-')
    combinedList = list(refTx = pre_report$ORF, testTx = queryTx)
    pre_report$txrevise_out =  indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
    
    pre_report$predictedStart = 'reconstructed'
    
    #pre_report = reconstructCDSstart(queryTx, knownCDS,
    #                                 refsequence,
    #                                 pre_report$txrevise_out,
    #                                 full.output = TRUE)
    output = utils::modifyList(output, list(predictedStart = pre_report$predictedStart))
    
    # return if CDS with new 5' do not contain a start codon
    if (is.na(pre_report$ORF[1])) {
      return(output)
    }
  } 
  
  # reconstruct CDS with insertion of alternative segments
  augmentedCDS = reconstructCDS(txrevise_out = pre_report$txrevise_out, 
                                fasta = refsequence, 
                                gene_id = gene_id, 
                                transcript_id = transcript_id)
  output = utils::modifyList(output, augmentedCDS)
  
  return(output)
}

