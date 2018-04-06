workflow1 <- function(knownCDS, queryTx, refseqeunce) {
  
  # prep output list
  output = list(transcript = queryTx,
             ORF = NA,
             AltSpliced_tx = NA,
             is_NMD = NA, 
             dist_to_lastEJ = NA, 
             annotatedStart = NA
             )
  
  # precheck for annotated start codon on query transcript and update output
  pre_report = testTXforStart(knownCDS, queryTx, full.output=TRUE)
  output = modifyList(output, pre_report[c("annotatedStart", "ORF")])
  
  # attempt to reconstruct CDS for transcripts with unannotated start
  if ((pre_report$annotatedStart == FALSE) |
      (pre_report$annotatedStart == TRUE & pre_report$firstexlen < 3)) {
    pre_report = reconstructCDSstart(knownCDS = knownCDS, 
                                     queryTx = queryTx,
                                     txrevise_out = pre_report$txrevise_out,
                                     refsequence = refseqeunce,
                                     full.output = TRUE)
    output = modifyList(output, pre_report["ORF"])

    # return if CDS with new 5' do not contain a start codon
    if (is.na(output$ORF[1])) {
      return(output)
    }
  } 

  # reconstruct CDS with insertion of alternative segments
  augmentedCDS = reconstructCDS(txrevise_out = pre_report$txrevise_out)
  output$AltSpliced_tx = ifelse(is.na(augmentedCDS[1]), FALSE, TRUE)
  
  # return if there is no alternative segments
  if (output$AltSpliced_tx == FALSE) {
    return(output)
  } else {
    NMDreport = testTXforNMD(augmentedCDS, refsequence, full.output = TRUE)
    output = modifyList(output, list(is_NMD = elementMetadata(NMDreport$stopcodons[1])$is_NMD, 
                                     dist_to_lastEJ = elementMetadata(NMDreport$stopcodons[1])$lastEJ_dist,
                                     ORF = NMDreport$ORF)
                        )
    return(output)
  }
}