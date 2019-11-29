

#' _ workflow: Test transcripts for NMD feature
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW.
#' This function will analyze the assembled transcripts for NMD features
#'
#' @param report_df 
#' dataframe of a list of assembled transcripts. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param basicTX_df 
#' dataframe of a list of reference CDS transcripts. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param inputExonsbyTx 
#' GRangesList of exon coordinates of assembled transcripts grouped by transcript name. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param basicExonsbyCDS 
#' GRangesList of exon coordinates of reference CDS grouped by transcript name. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param genome 
#' DNAstring containing genome sequence
#' @param PTC_dist 
#' Numerical value referring to minimium distance of premature stop codon to last exon junction
#' @param testNonClassicalNMD 
#' Boolean value (TRUE/FALSE) on whether to test for non-classical NMD features
#'
#' @return
#' updated form of report_df
#' @export
#'
#' @examples
runMain <- function(report_df, inputExonsbyTx, basicExonsbyCDS, 
                            basicExonsbyTx, genome, 
                            testforNMD = TRUE, PTC_dist = 50, 
                            testNonClassicalNMD = FALSE, 
                            testforAS = FALSE) {
  
  # this is an internal function for testing NMD features on report_df rowwise
  internalfunc = function(x) {
    
    # convert each line into a list so that elements 
    # can be referenced as thisline$_
    thisline = as.list(x)
    
    # Prepare GRanges and select best reference
    unpack[queryGRanges, basicTxGRanges, basicCDSGRanges, prepreport] = 
      getBestRef(thisline$Transcript_ID, thisline$Ref_transcript_ID,
                  thisline$Gene_ID, thisline$NMDer_ID,
                  inputExonsbyTx, basicExonsbyTx, basicExonsbyCDS) 
    
    # return if no overlap is found between query and any reference
    if(is.na(queryGRanges)){
      return(thisline)
    } else {
      thisline = utils::modifyList(thisline, prepreport)
    }
  
    # attempt to build ORF for query if absent
    if(thisline$ORF_found == FALSE){
      # attempt to build Open Reading Frame for query
      ORFreport = getORF(basicCDSGRanges, queryGRanges, genome)
      thisline = utils::modifyList(thisline, ORFreport)
    }

    # return if ORF is still not found 
    if(thisline$ORF_found == FALSE){
      # attempt to build Open Reading Frame for query
      return(thisline)
    }
    
    # if requested, test for NMD features and update line entry
    if (testforNMD == TRUE) {
      NMDreport = testNMD(thisline$ORF_considered, 
                          queryGRanges, 
                          PTC_dist, 
                          testNonClassicalNMD,
                          genome)
      thisline = utils::modifyList(thisline, NMDreport)
    }
    
    # if requested, classify alternative splicing events and update line entry
    if (testforAS == TRUE) {
      ASreport = classifyAS(queryGRanges, basicTxGRanges)
      
      
      
      if (testforNMD == FALSE) {
        ORF = NA
        is_NMD = NA
      } else {
        ORF = thisline$ORF_considered
        is_NMD = thisline$is_NMD
      }
      
      altevents = getASevents(basicTxGRanges, queryGRanges, 
                              testforNMD, ORF, is_NMD)
      
      thisline = utils::modifyList(thisline, altevents)
      
    }
    
    
    # update analyzed ORF coordinates into output
    if (thisline$ORF_found == TRUE) {
      thisline$ORF_considered = thisline$ORF_considered %>% as.data.frame()
    }
    rm(list = c('queryGRanges','basicCDSGRanges', 'basicTxGRanges'))
    return(thisline)

  }
  
  # run the above function on report_df
  report_df = report_df %>% 
    dplyr::rowwise() %>% 
    dplyr::do(data.frame(internalfunc(.), stringsAsFactors = FALSE)) 
  
  return(report_df)
}






