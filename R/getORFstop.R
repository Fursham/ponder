getORFstop <- function(query, fasta, fiveUTRlength){
  
  # prepare output list
  output = list(ORF_found = FALSE,
                threeUTRlength = 0)
  
  # append query GRanges to start from star codon, and retrieve seq
  queryCDS = resizeGRangesTranscripts(query, start = fiveUTRlength)
  queryseq = unlist(BSgenome::getSeq(fasta, queryCDS))
  
  # prepare a dict of stop codons for pattern matching
  list_stopcodons = Biostrings::DNAStringSet(c("TAA", "TAG", "TGA"))
  pdict_stopcodons = Biostrings::PDict(list_stopcodons)
  
  # search for in-frame stop codons
  stopcodons = Biostrings::matchPDict(pdict_stopcodons, queryseq) %>% 
    unlist() %>% 
    as.data.frame() %>%
    dplyr::filter(end %%3 ==0) %>%
    dplyr::arrange(start)
  
  # return if no in-frame stop codons are found
  if(nrow(stopcodons) == 0){
    return(output)
  } else{

    # retrieve 3UTR length and update output file
    firststopcodon = stopcodons[1,]
    threeUTRlength = length(queryseq) - firststopcodon$end
    output$threeUTRlength = threeUTRlength
    output$ORF_found = TRUE
    
    return(output)
  }
  
}