getORF <- function(query, CDS, fasta) {
  
  # prepare output list
  output = list(ORF_considered = as.character(NA),
                ORF_start = as.character('Not found'),
                fiveUTRlength = 0,
                threeUTRlength = 0,
                ORF_found = FALSE)
  
  queryTx = query[[1]]
  knownCDS = CDS[[1]]
  mcols(queryTx)$transcript_id = names(query)
  # attempt to find an aligned start codon
  report = getORFstart(queryTx, knownCDS, fasta)
  output = utils::modifyList(output, report) #update output
  
  # return if no start codon is found
  if(output$ORF_start == 'Not found'){
    return()
  }
  
  # attempt to search for an in-frame stop codon
  report = getORFstop(queryTx,  fasta, output$fiveUTRlength)
  output = utils::modifyList(output, report) #update output
  
  # return if no stop codon is found
  if(output$ORF_found == FALSE){
    return()
  }
  
  # build new ORF Granges
  report = getORFranges(queryTx, output$fiveUTRlength, output$threeUTRlength)
  output = utils::modifyList(output, report) #update output
  #return(output[c('ORF_considered', 'ORF_start', 'ORF_found')])
  return(output$ORF_considered)
  
}