#' This function will test if the CDS in query is an NMD substrate
#'
#' @param queryCDS is a GRanges Object
#' @param refsequence is a sequence containing TwoBitFile
#'
#' @return Granges
#' @export
testNMDtx <- function(queryCDS, refsequence){
  
  queryStrand = as.character(strand(queryCDS))[1]
  
  # prepares query seq and exon boundaries depending on strand
  if (queryStrand == '-') {
    queryCDS = sort(queryCDS)
    thisqueryseq = getSeq(mmus_dna, queryCDS)
    thisqueryseq = reverseComplement(unlist(thisqueryseq))
    exon_boundaries = cumsum(rev(width(queryCDS)))
  } else if(queryStrand == '+') {
    thisqueryseq = getSeq(mmus_dna, queryCDS)
    exon_boundaries = cumsum(width(queryCDS))
  }
  
  # prepare a dict of stop codons for pattern matching
  stopcodons = DNAStringSet(c("TAA", "TAG", "TGA"))
  pdict_stopcodons = PDict(stopcodons)
  
  # search for in-frame stop codons
  matches = matchPDict(pdict_stopcodons, thisqueryseq)
  combinedmatches = do.call("c", matches)
  stop_codon = combinedmatches[end(combinedmatches) %% 3 == 0,]
  PTC = stop_codon[end(stop_codon) != sum(width(queryCDS))]
  return(PTC)
  PTC_rep <- ifelse(!is.null(PTC), 'Yes', 'No')  # need to improve on this
  
  # calculate distance of stop codon to last exon junction
  lastEJ = head(tail(exon_boundaries, n=2), n=1)
  dist_stop_to_lastEJ = end(stop_codon) - lastEJ
  
  metadata = dplyr::data_frame(lastEJ_dist = dist_stop_to_lastEJ) %>%
    dplyr::mutate(is_NMD = ifelse(lastEJ_dist < -50, TRUE, FALSE)) %>%
    as.data.frame()
  elementMetadata(stop_codon) = metadata
  return(stop_codon)
  
  # outputs (1) whether tx contain PTC and (2) distance of (all) stops to last EJ
  # more outputs would be desirable etc exon structure considered
  output = list(PTC_rep, dist_stop_to_lastEJ)
  return(output)
}
