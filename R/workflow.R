workflow1 <- function(knownCDS, queryTx, refseqeunce) {
  
  # precheck for annotated start codon on query transcript
  pre_report = testTXforStart(knownCDS, queryTx)
  
  # attempt to reconstruct CDS for transcripts with unannotated start
  if (pre_report$annotatedStart == FALSE) {
    pre_report = reconstructCDSstart(knownCDS = knownCDS, 
                                     queryTx = queryTx,
                                     txrevise_out = pre_report$txrevise_out,
                                     refsequence = refseqeunce)
    queryTx = pre_report$newCDS
    
    # return if CDS with new 5' do not contain a start codon
    if (is.na(pre_report$txrevise_out[1])) {
      return(list(transcript = queryTx,
                  AltSpliced_tx = FALSE,
                  is_NMD = NA, 
                  dist_to_lastEJ = NA, 
                  annotatedStart = pre_report$annotatedStart))
    }
  } 

  # reconstruct CDS with insertion of alternative segments
  augmentedCDS = reconstructCDS(txrevise_out = pre_report$txrevise_out, 
                                annotatedStart = pre_report$annotatedStart)
  
  # return if there is no
  if (is.na(augmentedCDS[[1]][1])) {
    return(list(transcript = queryTx,
                AltSpliced_tx = FALSE,
                is_NMD = NA, 
                dist_to_lastEJ = NA, 
                annotatedStart = pre_report$annotatedStart)
           )
  }
  
  NMDreport = testTXforNMD(augmentedCDS[[1]], refsequence)
  return(list(transcript = augmentedCDS[[1]], 
              AltSpliced_tx = TRUE,
              is_NMD = elementMetadata(NMDreport[1])$is_NMD, 
              dist_to_lastEJ = as.numeric(elementMetadata(NMDreport[1])$lastEJ_dist), 
              annotatedStart = pre_report$annotatedStart)
         )
  
  
}