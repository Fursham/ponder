getORFstop <- function(query, fasta, fiveUTRlength){
  
  # prepare output list
  output = list(ORF_found = FALSE,
                threeUTRlength = 0)
  
  # append query GRanges to start from star codon, and retrieve seq
  queryCDS = resizeTranscripts(query, headlength = fiveUTRlength)
  queryseq = unlist(Biostrings::getSeq(fasta, queryCDS))
  
  # prepare a dict of stop codons for pattern matching
  list_stopcodons = Biostrings::DNAStringSet(c("TAA", "TAG", "TGA"))
  pdict_stopcodons = Biostrings::PDict(list_stopcodons)
  
  # search for in-frame stop codons
  stopcodons = Biostrings::matchPDict(pdict_stopcodons, queryseq) %>% 
    unlist()
  stopcodons = sort(stopcodons[end(stopcodons) %% 3 == 0,])
  
  # return if no in-frame stop codons are found
  if(length(stopcodons) == 0){
    return(output)
  } else{

    # retrieve 3UTR length and update output file
    firststopcodon = stopcodons[1]
    threeUTRlength = length(queryseq) - end(firststopcodon)
    output$threeUTRlength = threeUTRlength
    output$ORF_found = TRUE
    
    return(output)
  }
  
}