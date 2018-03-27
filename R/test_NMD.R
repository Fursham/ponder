#' Generate augmented CDSs of spliced transcripts
#'
#' @param knownCDS GRanges object containing reference CDS for the gene in query
#' @param queryTx GRanges object containing exon structure of query transcript
#'
#' @return GRangesList object containing a list of new CDS structures or 
#' NULL if no unique internal exons are found
#' @author Fursham Hamid
#' @export
#'
augmentCDS <- function(knownCDS, queryTx){

  # combine both GRanges objects into a list and identify segments which are different
  combinedList = list(refTx = knownCDS, testTx = queryTx)
  diffSegments = indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
  
  if (length(diffSegments$refTx$contained[diffSegments$refTx$contained == TRUE]) == 0 & 
      length(diffSegments$testTx$contained[diffSegments$testTx$contained == TRUE]) == 0) {
    return(NULL)
  }
  # remove internal segments in CDS which are absent in test transcripts,
  # add internal segments to CDS which are present in test transcripts
  ## future development, to also test first and last exons
  augmentedTx = knownCDS
  if (!is.null(diffSegments$refTx$contained)){
    augmentedTx = knownCDS[knownCDS != diffSegments$refTx[diffSegments$refTx$contained == TRUE]]
    augmentedTx = reduce(augmentedTx)
  }
  if (!is.null(diffSegments$testTx$contained)){
    augmentedTx = sort(unlist(append(
      reduce(augmentedTx), 
      reduce(diffSegments$testTx[diffSegments$testTx$contained == TRUE]))),
      decreasing = as.character(strand(knownCDS))[1] == '-')
  }

  return(augmentedTx) 
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
    thisqueryseq = Biostrings::getSeq(mmus_dna, queryCDS) %>%
      unlist(.)
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
