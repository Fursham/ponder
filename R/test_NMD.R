#' Generate augmented CDSs of spliced transcripts
#'
#' @param knownCDS GRanges object containing reference ORF for the gene in query
#' @param queryTx GRanges object containing exon structure of query transcript
#'
#' @return GRangesList object containing a list of new CDS structures
#' @export
#'
augmentCDS <- function(knownCDS, queryTx){

  # merge both GRanges transcripts
  combinedList = removeMetadata(list(refTx = knownCDS, testTx = queryTx))
  
  
  # identify segments which are different
  diffSegments = indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
  
  # we need to find a way to remove segments found in ref and add segments found in test
  augmentedTx = setdiff(knownCDS, diff$ENSMUST00000029780[2,])
  
  
  return() 
}




#' Test query CDS for NMD features
#'
#' @param queryCDS GRanges Object containing CDS structure
#' @param refsequence sequence file build from AnnotationHub
#' 
#' @return A GRanges object containing a list of in-frame stop codons, its
#' coordinate, its distance to last EJC and prediction of its NMD nature
#' @import Biostrings 
#' @author Fursham Hamid
#' @export
testTxforNMD <- function(queryCDS, refsequence){
  
  queryStrand = as.character(strand(queryCDS))[1]
  
  # prepares query seq and exon boundaries depending on strand
  if (queryStrand == '-') {
    queryCDS = sort(queryCDS)
    thisqueryseq = Biostrings::getSeq(mmus_dna, queryCDS)
    thisqueryseq = Biostrings::reverseComplement(unlist(thisqueryseq))
    exon_boundaries = cumsum(rev(width(queryCDS)))
  } else if(queryStrand == '+') {
    thisqueryseq = Biostrings::getSeq(mmus_dna, queryCDS)
    exon_boundaries = cumsum(width(queryCDS))
    
  }

  # prepare a dict of stop codons for pattern matching
  list_stopcodons = Biostrings::DNAStringSet(c("TAA", "TAG", "TGA"))
  pdict_stopcodons = Biostrings::PDict(list_stopcodons)
  
  # search for in-frame stop codons
  allmatches = Biostrings::matchPDict(pdict_stopcodons, thisqueryseq) # check whether this list is named
  combinedmatches = unlist(allmatches)
  inframe_stopcodons = combinedmatches[end(combinedmatches) %% 3 == 0,]

  
  ## this line is no longer necessary
  # PTC = inframe_stopcodons[end(inframe_stopcodons) != sum(width(queryCDS))]
  
  # calculate distance of stop codon to last exon junction
  lastEJ = head(tail(exon_boundaries, n=2), n=1)
  dist_stop_to_lastEJ = end(inframe_stopcodons) - lastEJ
  
  # update GRanges object 'inframe_stopcodons' with more information
  metadata = dplyr::data_frame(lastEJ_dist = dist_stop_to_lastEJ) %>%
    dplyr::mutate(is_NMD = ifelse(lastEJ_dist < -50, TRUE, FALSE)) %>%
    as.data.frame()
  elementMetadata(inframe_stopcodons) = metadata
  return(inframe_stopcodons)
}
