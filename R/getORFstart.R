getORFstart <- function(query, refCDS, fasta){
  
  output = list(ORF_start = as.character('Not found'),
                fiveUTRlength = 0)
  
  # get coord of start codon on ref
  startcodon = resizeTranscripts(refCDS, end = sum(width(refCDS))-3)
  strand = as.character(strand(query))[1]
  
  # annotated start is within query
  if(startcodon %within% query & length(startcodon) == 1){
    disjoint = BiocGenerics::append(query,startcodon) %>%
      GenomicRanges::disjoin(with.revmap = T) %>%
      sort(decreasing = strand == '-')
    mcols(disjoint)$cumsum = cumsum(width(disjoint))
      
    startcodonindex = min(which(lengths(mcols(disjoint)$revmap) == 2)) -1
    fiveUTRlength = mcols(disjoint)$cumsum[startcodonindex]
    
    output$ORF_start = 'Annotated'
    output$fiveUTRlength = fiveUTRlength
    
    return(output)
  } 
  # attempt to find in-frame ATG
  else {
    refsequence = unlist(Biostrings::getSeq(fasta, refCDS))
    startcodons = matchPattern('ATG', refsequence) %>% ranges()
    inframestarts = startcodons[end(startcodons) %% 3 == 0 & 
                                  start(startcodons) != 1]
    
    if(length(inframestarts) == 0){
      return(output)
    } else {

      inframestartsingranges = do.call('c', base::mapply(function(x,y){
        start = x - 1
        end = length(refsequence) - y
        startcodoninGRanges = resizeTranscripts(refCDS, start, end)
        return(startcodoninGRanges)
      }, start(inframestarts), end(inframestarts)))
      
      if(any(inframestartsingranges %within% query)){
        startgranges = inframestartsingranges[inframestartsingranges %within% query][1]
        disjoint = BiocGenerics::append(query,startgranges) %>%
          GenomicRanges::disjoin(with.revmap = T) %>%
          sort(decreasing = strand == '-')
        mcols(disjoint)$cumsum = cumsum(width(disjoint))
        
        startcodonindex = match(2,lengths(mcols(disjoint)$revmap)) -1
        fiveUTRlength = mcols(disjoint)$cumsum[startcodonindex]
        
        output$ORF_start = 'Predicted'
        output$fiveUTRlength = fiveUTRlength
        
        return(output)
      
      
      } else{
        return(output)
      }
    }
  }
}
