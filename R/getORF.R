getORF <- function(knownCDS, queryTx, fasta) {
  
  # prepare output list
  output = list(ORF_considered = as.character(NA),
                ORF_start = as.character('Not found'),
                fiveUTRlength = 0,
                threeUTRlength = 0,
                ORF_found = FALSE)
  
  # attempt to find an aligned start codon
  report = getORFstart(queryTx, knownCDS, fasta)
  output = utils::modifyList(output, report) #update output
  
  # return if no start codon is found
  if(output$ORF_start == 'Not found'){
    return(output[c('ORF_considered', 'ORF_start', 'ORF_found')])
  }
  
  # attempt to search for an in-frame stop codon
  report = getORFstop(queryTx,  fasta, output$fiveUTRlength)
  output = utils::modifyList(output, report) #update output
  
  # return if no stop codon is found
  if(output$ORF_found == FALSE){
    return(output[c('ORF_considered', 'ORF_start', 'ORF_found')])
  }
  
  # build new ORF Granges
  report = getCDSranges(queryTx, output$fiveUTRlength, output$threeUTRlength)
  output = utils::modifyList(output, report) #update output
  return(output[c('ORF_considered', 'ORF_start', 'ORF_found')])
  
}