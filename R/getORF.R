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
#' @param gene_id 
#' @param transcript_id 
#'
#' @return
#'
#' @examples
getORF <- function(knownCDS, queryTx, refsequence, gene_id, transcript_id) {
  
  # prep output list
  output = list(ORF_considered = as.character(NA),
                ORF_start = as.character('Not found'),
                ORF_found = FALSE)
  
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
    
    pre_report = reconstructCDSstart(queryTx, knownCDS,
                                     refsequence,
                                     pre_report$txrevise_out,
                                     full.output = TRUE)
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

