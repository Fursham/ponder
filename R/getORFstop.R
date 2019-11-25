getORFstop <- function(query, fasta, fiveUTRlength){
  
  output = list(ORF_found = FALSE,
                threeUTRlength = 0)
  
  strand = as.character(strand(query))[1]
  queryCDS = resizeTranscripts(query, start = fiveUTRlength)
  queryseq = unlist(Biostrings::getSeq(fasta, queryCDS))
  
  # prepare a dict of stop codons for pattern matching
  list_stopcodons = Biostrings::DNAStringSet(c("TAA", "TAG", "TGA"))
  pdict_stopcodons = Biostrings::PDict(list_stopcodons)
  
  # search for in-frame stop codons
  stopcodons = Biostrings::matchPDict(pdict_stopcodons, queryseq) %>% 
    unlist()
  stopcodons = sort(stopcodons[end(stopcodons) %% 3 == 0,])
  
  if(length(stopcodons) == 0){
    return(output)
  } else{

    firststopcodon = stopcodons[1]
    threeUTRlength = length(queryseq) - end(firststopcodon)
    output$threeUTRlength = threeUTRlength
    output$ORF_found = TRUE
    
    return(output)
  }
  
}