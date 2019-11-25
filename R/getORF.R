getORF <- function(knownCDS, queryTx, fasta) {
  
  # prep output list
  output = list(ORF_considered = as.character(NA),
                ORF_start = as.character('Not found'),
                fiveUTRlength = 0,
                threeUTRlength = 0,
                ORF_found = FALSE)
  
  # attempt to get start codon
  report = getORFstart(queryTx, knownCDS, fasta)
  output = utils::modifyList(output, report)
  
  # return if no start codon is found
  if(output$ORF_start == 'Not_found'){
    return(output[c('ORF_considered', 'ORF_start', 'ORF_found')])
  }
  
  # attempt to get stop codon
  report = getORFstop(queryTx,  fasta, output$fiveUTRlength)
  output = utils::modifyList(output, report)
  
  # return if no stop codon is found
  if(output$ORF_found == FALSE){
    return(output[c('ORF_considered', 'ORF_start', 'ORF_found')])
  }
  
  report = getCDSranges(queryTx, output$fiveUTRlength, output$threeUTRlength)
  output = utils::modifyList(output, report)
  return(output[c('ORF_considered', 'ORF_start', 'ORF_found')])
  
}