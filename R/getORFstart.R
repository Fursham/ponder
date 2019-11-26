getORFstart <- function(query, refCDS, fasta){
  
  # prepare output list
  output = list(ORF_start = as.character('Not found'),
                fiveUTRlength = 0)
  
  # get coord of start codon on reference and strand info
  startcodon = resizeGRangesTranscripts(refCDS, end = sum(width(refCDS))-3)
  strand = as.character(strand(query))[1]
  
  # if query containg annotated start codon:
  if(startcodon %within% query & length(startcodon) == 1){
    
    # this part of the code will calculate the length of 5'UTR
    #   disjoin will break query GRanges into sub GRanges, based
    #   on the position of the start codon.
    #   we can then find length of the upstream/downstream segments 
    #   of the break
    disjoint = BiocGenerics::append(query,startcodon) %>%
      GenomicRanges::disjoin(with.revmap = T) %>%
      sort(decreasing = strand == '-')
    
    # calculates cumulative sum of segment
    cumsumwidth = cumsum(width(disjoint))
      
    # retrieve index of segment upstream of start codon and return its cumsumwidth
    startcodonindex = min(which(lengths(mcols(disjoint)$revmap) == 2)) -1
    fiveUTRlength = cumsumwidth[startcodonindex]
    
    # update output list
    output$ORF_start = 'Annotated'
    output$fiveUTRlength = fiveUTRlength
    
    return(output)
  } 
  # if annotated start is not found, attempt to find upstream-most internal ATG
  else {
    
    # get sequence of ref, find all internal inframe-ATG
    refsequence = unlist(Biostrings::getSeq(fasta, refCDS)) # 
    startcodons = Biostrings::matchPattern('ATG', refsequence) %>% ranges()
    inframestarts = startcodons[end(startcodons) %% 3 == 0 & 
                                  start(startcodons) != 1]
    
    # return if no internal ATG is found
    if(length(inframestarts) == 0){
      return(output)
    } else {

      # This function attempts to map the XStringViews output back to refGRanges
      inframestartsingranges = do.call('c', base::mapply(function(x,y){
        start = x - 1
        end = length(refsequence) - y
        startcodoninGRanges = resizeGRangesTranscripts(refCDS, start, end)
        return(startcodoninGRanges)
      }, start(inframestarts), end(inframestarts)))
      
      # obtain 5'UTR length if query contain any of the inframe ATG
      if(any(inframestartsingranges %within% query)){
        
        firststartgranges = inframestartsingranges[inframestartsingranges %within% query][1]
        disjoint = BiocGenerics::append(query,firststartgranges) %>%
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
