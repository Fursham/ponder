#' Test transcripts for annotated start codons
#'
#' @description This function will test whether the query transcript contain an annotated 
#' Start codon from a reference CCDS of the gene
#' 
#' @param knownCDS A GRanges object containing CCDS information of a gene family
#' @param queryTx A GRanges object containing exon structure of a query transcript from
#' the same gene family
#'
#' @return 
#' A list of GRanges objects similar to 'indentifyAddedRemovedRegions()' output. 
#' Fourth object in the list flags the transcript for annotated start codon
#' ?indentifyAddedRemovedRegions for more information
#' @export
#' 
#' @author Fursham Hamid
#' 
#' @examples
#' 
#' testTXforStart(ptbp2_testData$noNMD, ptbp2_testData$NMD)
#' 
#' 
testTXforStart <- function(knownCDS, queryTX) {
  
  # combine both GRanges objects into a list and identify segments which are different
  combinedList = list(refTx = knownCDS, testTx = queryTX)
  diffSegments = indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
  
  # test if the reference CDS contain any unique upstream segments
    # if test return FALSE, it means that the test transcript bear the annotated Start
    # if test return TRUE, it means that the test transcript bear a different 5' structure
  if (length(diffSegments$refTx$upstream[diffSegments$refTx$upstream == TRUE]) == 0) {
    output = list(txrevise_out = diffSegments, annotatedStart = TRUE)
  } else {
    output = list(txrevise_out = diffSegments, annotatedStart = FALSE)
  }
  return(output)
}


#' Reconstruct CDS with alternate 5' exons
#'
#' @description This function will reconstruct a CDS sequence with a different 5' end 
#' taken from a query transcript from the same gene
#' 
#' @param knownCDS A GRanges object containing CCDS information of a gene family
#' @param queryTx A GRanges object containing exon structure of a query transcript from
#' the same gene family
#' @param txrevise_out GRangesobject from 'indentifyAddedRemovedRegions()' 
#' or 'testTXforStart()' output. ?indentifyAddedRemovedRegions or ?testTXforStart
#' for more information
#' @param refsequence sequence file build from AnnotationHub
#' 
#' @return 
#' A list containing (1) a GRanges object with a reconstructed CDS, 
#' (2) GRangesobject from 'indentifyAddedRemovedRegions()' and 
#' (3) a flag for non-annotated Start codon.
#' If reconstructed CDS do not contain an in-frame Start codon with 3' end of the coding sequence,
#' function will return NA for (2)
#' 
#' @export
#' 
#' @author Fursham Hamid
#'
#' @examples
#' 
#' reconstructCDSstart(ptbp2_testData$noNMD, ptbp2_testData$diffstart, refsequence = mmus_dna)
#' 
reconstructCDSstart <- function(knownCDS, queryTx, txrevise_out, refsequence) {
  
  # check if txrevise_out input is provided, and build one if not.
  if (missing(knownCDS) | missing(queryTx)) {
    stop('Please provide input GRanges objects')
  } else if (missing(txrevise_out)) {
    combinedList = list(refTx = knownCDS, testTx = queryTx)
    diffSegments = indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
  } else {
    diffSegments = txrevise_out
  }
  
  # obtain strand information and anchor exon
  queryStrand = as.character(strand(diffSegments[[3]]))[1]
  anchor_exon = diffSegments[[3]][1]
  
  # choose downstream exon as anchor if first shared exon overlaps with first CDS exons
  if ((is.na(as.matrix(findOverlaps(knownCDS[1], diffSegments[[3]][1]))[1]) == TRUE)) {
    anchor_exon = diffSegments[[3]][2]
  } 
 
  # reconstruct CDS with appended new 5' end
  CDSsegment_to_be_retained = knownCDS[which(knownCDS %in% anchor_exon):length(knownCDS)]
  TxSegment_to_be_retained = queryTx[1:(which(queryTx %in% anchor_exon)-1)]
  
  reconstructedTx = sort(unlist(append(
    reduce(TxSegment_to_be_retained), 
    reduce(CDSsegment_to_be_retained))),
    decreasing =  queryStrand == '-')
  
  # prepares query seq and exon boundaries depending on strand
  if (queryStrand == '-') {
    queryCDS = sort(reconstructedTx)
    thisqueryseq = Biostrings::getSeq(mmus_dna, queryCDS)
    thisqueryseq = Biostrings::reverseComplement(unlist(thisqueryseq))
  } else if(queryStrand == '+') {
    thisqueryseq = Biostrings::getSeq(mmus_dna, reconstructedTx) %>%
      unlist(.)
  }
  

  # find all in-frame start codons in reconstructed transcript
  StartCodons = Biostrings::matchPattern(Biostrings::DNAString("ATG"), thisqueryseq)
  inFrameStartCodons = StartCodons[(length(reconstructedTx) - (end(StartCodons))) %%3 == 0,]
  
  if (length(inFrameStartCodons) > 0) {

    # obtain txrevise output
    combinedList = list(refTx = reconstructedTx, testTx = queryTx)
    diffSegments = indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
    
    out = list(newCDS = reconstructedTx, txrevise_out = diffSegments, annotatedStart = FALSE)
  } else {
    out = list(newCDS = reconstructedTx, txrevise_out = NA, annotatedStart = FALSE)
  }
  return(out)
}


#' Reconstruct ORF with alternative segments
#'
#' @description 
#' This function will generate a new ORF of a gene family with 
#' added segments from a query transcript. INPUT...
#' 
#' @param knownCDS A GRanges object containing CCDS information of a gene family
#' @param queryTx A GRanges object containing exon structure of a query transcript from
#' the same gene family
#' @param txrevise_out 
#' Optional. 
#' GRangesobject from 'indentifyAddedRemovedRegions()' 
#' or 'testTXforStart()' output. ?indentifyAddedRemovedRegions or ?testTXforStart
#' for more information.
#' Arguments knownCDS and queryTx is not mandatory if this argument is provided
#' @param annotatedStart A 
#' 
#' 
#' @return
#' A GRanges Object with new ORF information if transcript contain unique segments, or
#' NA if transcript is identical to reference CDS
#' @export
#' 
#' @author Fursham Hamid
#'
#' @examples
#' 
#' reconstructCDS(_)
#' 
#' 
#' 
#' 
#' 
reconstructCDS <- function(knownCDS, queryTx, txrevise_out, annotatedStart = TRUE){

  # check if txrevise_out input is provided, and build one if not.
  if (missing(txrevise_out)) {
    if (missing(knownCDS) | missing(queryTx)) {
      stop('Please provide input GRanges objects')
    }
    combinedList = list(refTx = knownCDS, testTx = queryTx)
    diffSegments = indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
  } else {
    diffSegments = txrevise_out
  }
  
  # return NULL if there is no unique internal and downstream segments
    # special case
  if (length(diffSegments[[1]]$contained[diffSegments[[1]]$contained == TRUE]) == 0 &
      length(diffSegments[[1]]$downstream[diffSegments[[1]]$downstream == TRUE]) == 0 &
      length(diffSegments[[2]]$contained[diffSegments[[2]]$contained == TRUE]) == 0
      ) {
    if(annotatedStart == TRUE){
      output = list(NA, annotatedStart = TRUE)
    } else {
      output = list(NA, annotatedStart = FALSE)
    }
    return(output)
  } 

  augmentedTx = sort(reduce(unlist(append(
    diffSegments$shared_exons, 
    reduce(diffSegments[[2]][diffSegments[[2]]$upstream != TRUE])))),
    decreasing = as.character(strand(diffSegments$shared_exons))[1] == '-')

  return(list(augmentedTx))
}


#' Test ORF for NMD
#' 
#' @description 
#' This function will check for NMD features for an input transcript
#'
#' @param queryCDS GRanges Object containing CDS structure
#' @param refsequence sequence file build from AnnotationHub
#' 
#' @return A GRanges object containing a list of NMD-inducible stop codons
#' 
#' @import Biostrings
#' @author Fursham Hamid
#' @export
#' 
#' @examples
#'
#' # Return a GRanges Object with 1 stop codon
#' testTxforNMD(ptbp2_testData$noNMD, mmus_dna)
#' 
#' # Return a GRanges Object with multiple premature stop codons, 
#' # 5 of which is >50 from lastEJC
#' testTxforNMD(ptbp2_testData$NMD, mmus_dna)
#' 
testTXforNMD <- function(queryCDS, refsequence){
  
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

  # check first codon for ATG and append to nearest ATG if it is not
  if (thisqueryseq[1:3] != Biostrings::DNAString("ATG")) {
    StartCodons = Biostrings::matchPattern(Biostrings::DNAString("ATG"), thisqueryseq)
    nearestStart = start(StartCodons[1])
    thisqueryseq = thisqueryseq[nearestStart:length(thisqueryseq)]
    
    UTRlength = nearestStart - 1
  } else {
    UTRlength = 0
  }
  # check for at least a start and end is an ORF
  
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
  lastEJ = head(tail(exon_boundaries, n=2), n=1) - UTRlength
  dist_stop_to_lastEJ = end(inframe_stopcodons) - lastEJ
  
  # update GRanges object 'inframe_stopcodons' with more information
  metadata = dplyr::data_frame(lastEJ_dist = dist_stop_to_lastEJ) %>%
    dplyr::mutate(is_NMD = ifelse(lastEJ_dist < -50, TRUE, FALSE)) %>%
    as.data.frame()
  elementMetadata(inframe_stopcodons) = metadata
  inframe_stopcodons = inframe_stopcodons[order(elementMetadata(inframe_stopcodons)$lastEJ_dist)]
  
  #inframe_stopcodons = inframe_stopcodons[elementMetadata(inframe_stopcodons)$is_NMD == TRUE] %>%
  #  .[order(elementMetadata(.)$lastEJ_dist)]
  return(inframe_stopcodons)
}
