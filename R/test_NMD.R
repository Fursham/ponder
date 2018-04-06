#' Test transcripts for annotated start codons
#'
#' @description This function will test whether the query transcript contain an annotated 
#' Start codon from a reference CCDS of the gene
#' 
#' @usage 
#' testTXforStart(knownCDS, queryTX, full.output = FALSE)
#' 
#' @param knownCDS A GRanges object containing CCDS information of a gene family
#' @param queryTx A GRanges object containing exon structure of a query transcript from
#' the same gene family
#' @param full.output If TRUE, function will output other information 

#' @return
#' By default, function will return a TRUE/FALSE boolean (annotatedStart)
#' 
#' If full.output is TRUE, function will output other information which include:
#' (1) the length of the first coding exon (firstexlen), 
#' (2) queryTx GRanges with 5' end appended to the start codon  (ORF),  
#' (3) a list of GRanges objects similar to 'indentifyAddedRemovedRegions()' output (txrevise_out). 
#' ?indentifyAddedRemovedRegions for more information
#' 
#' @export
#' 
#' @author Fursham Hamid
#' 
#' @examples
#' 
#' testTXforStart(ptbp2_testData$noNMD, ptbp2_testData$NMD)
#' testTXforStart(ptbp2_testData$noNMD, ptbp2_testData$diffstart)
#' 
#' 
testTXforStart <- function(knownCDS, queryTX, full.output = FALSE) {
  
  # prepare output of function
  output = list(annotatedStart = NA,
                txrevise_out = NA,
                ORF = NA,
                firstexlen = NA)

  # combine both GRanges objects into a list and identify segments which are different
  combinedList = list(refTx = knownCDS, testTx = queryTX)
  diffSegments = indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
  
  # test if the reference CDS contain any unique upstream segments
    # if test return FALSE, it means that the test transcript bear the annotated Start
    # if test return TRUE, it means that the test transcript bear a different 5' structure
  
  if (length(diffSegments$refTx$upstream[diffSegments$refTx$upstream == TRUE]) == 0) {
    
    lenfirstsharedexon = width(diffSegments$shared_exons[1])
    upUTRsize = sum(width(diffSegments$testTx[diffSegments$testTx$upstream == TRUE]))
    setORF = resizeTranscripts(queryTX, start = upUTRsize)
    output = modifyList(output, 
                        list(txrevise_out = diffSegments, 
                             ORF = setORF, 
                             annotatedStart = TRUE, 
                             firstexlen = lenfirstsharedexon))
    } else {
      output = modifyList(output, 
                          list(txrevise_out = diffSegments, 
                               annotatedStart = FALSE))
  }
  ifelse(full.output == TRUE, return(output), return(output["annotatedStart"]))
}


#' Reconstruct CDS with alternate 5' exons
#'
#' @description This function will reconstruct a CDS sequence with a different 5' end 
#' taken from a query transcript from the same gene
#' 
#' @usage reconstructCDSstart(knownCDS, queryTx, txrevise_out, refsequence, full.output = FALSE)
#' 
#' @param knownCDS A GRanges object containing CCDS information of a gene family
#' @param queryTx A GRanges object containing exon structure of a query transcript from
#' the same gene family
#' @param txrevise_out GRangesobject from 'indentifyAddedRemovedRegions()' 
#' or 'testTXforStart()' output. ?indentifyAddedRemovedRegions or ?testTXforStart
#' for more information
#' @param refsequence sequence file build from AnnotationHub
#' @param full.output if TRUE, function will output additional information
#' 
#' @return 
#' Default output (ORF),
#' If reconstructed CDS contain an in-frame ATG, function will output a GRanges Object describing the
#' reconstructed CDS. Else, function will output 'NA'
#' 
#' If full.output == TRUE,
#' function will also return a list of GRanges objects similar to 'indentifyAddedRemovedRegions()' output (txrevise_out). 
#' ?indentifyAddedRemovedRegions for more information
#' 
#' @export
#' 
#' @author Fursham Hamid
#'
#' @examples
#' 
#' reconstructCDSstart(ptbp2_testData$noNMD, ptbp2_testData$diffstart, refsequence = mmus_dna)
#' 
reconstructCDSstart <- function(knownCDS, queryTx, txrevise_out, refsequence, full.output = FALSE) {
  
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
    
    # append 5' end of reconstructed transcript
    upUTRsize = start(inFrameStartCodons[1]) - 1
    setORF = resizeTranscripts(reconstructedTx, start = upUTRsize)
    
    # obtain txrevise output
    combinedList = list(refTx = setORF, testTx = queryTx)
    diffSegments = indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
    out = list(ORF = setORF, txrevise_out = diffSegments)
  } else {
    out = list(ORF = NA, txrevise_out = NA)
  }
  ifelse(full.output == TRUE, return(out), return(out["ORF"]))
}


#' Reconstruct ORF with alternative internal and downstream segments
#'
#' @description 
#' This function will generate a new ORF of a gene family with 
#' added segments from a query transcript. 
#' 
#' Do note that this function, if ran as a stand-alone, 
#' will not insert/remove first exons. 
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
#' reconstructCDS(ptbp2_testData$noNMD, ptbp2_testData$exons$ENSMUST00000197833)
#' 
#' 
#' 
#' 
#' 
reconstructCDS <- function(knownCDS, queryTx, txrevise_out){

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
    return(NA)
  } 

  augmentedTx = sort(reduce(unlist(append(
    diffSegments$shared_exons, 
    reduce(diffSegments[[2]][diffSegments[[2]]$upstream != TRUE])))),
    decreasing = as.character(strand(diffSegments$shared_exons))[1] == '-')

  return(augmentedTx)
}


#' Test ORF for NMD
#' 
#' @description 
#' This function will check for NMD features for an input ORF
#' 
#' Do not that input query has to begin with an ATG
#'
#' @param queryCDS GRanges Object containing CDS structure
#' @param refsequence sequence file build from AnnotationHub
#' @param full.output If TRUE, function will output other information 
#' 
#' @return 
#' A list containing: 
#' (1) GRanges object containing a list of in-frame stop codons and its metadata
#' 
#' if full.output == TRUE, function will also return
#' (2) An ORF which ends at the most upstream in-frame stopcodon
#' 
#' @import Biostrings
#' @author Fursham Hamid
#' @export
#' 
#' @examples
#' 
#' testTXforNMD(ptbp2_testData$noNMD, mmus_dna)
#' testTXforNMD(ptbp2_testData$NMD, mmus_dna)
#' 
testTXforNMD <- function(queryCDS, refsequence, minmPTCdist = -50, full.output = FALSE){
  
  queryStrand = as.character(strand(queryCDS))[1]
  exon_boundaries = cumsum(width(queryCDS))
  
  # prepares query seq and exon boundaries depending on strand
  if (queryStrand == '-') {
    setqueryCDS = sort(queryCDS)
    thisqueryseq = Biostrings::getSeq(mmus_dna, setqueryCDS)
    thisqueryseq = Biostrings::reverseComplement(unlist(thisqueryseq))
  } else if(queryStrand == '+') {
    thisqueryseq = Biostrings::getSeq(mmus_dna, queryCDS) %>%
      unlist(.)
  }
  
  # check first codon for ATG and append to nearest ATG if it is not
  #if (thisqueryseq[1:3] != Biostrings::DNAString("ATG")) {
  #  StartCodons = Biostrings::matchPattern(Biostrings::DNAString("ATG"), thisqueryseq)
  #  nearestStart = start(StartCodons[1])
  #  thisqueryseq = thisqueryseq[nearestStart:length(thisqueryseq)]
    
  #  UTRlength = nearestStart - 1
  #} else {
  #  UTRlength = 0
  #}
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
  lastEJ = head(tail(exon_boundaries, n=2), n=1)
  dist_stop_to_lastEJ = end(inframe_stopcodons) - lastEJ
  
  # update GRanges object 'inframe_stopcodons' with more information
  metadata = dplyr::data_frame(lastEJ_dist = dist_stop_to_lastEJ) %>%
    dplyr::mutate(is_NMD = ifelse(lastEJ_dist < minmPTCdist, TRUE, FALSE)) %>%
    as.data.frame()
  elementMetadata(inframe_stopcodons) = metadata
  inframe_stopcodons = inframe_stopcodons[order(elementMetadata(inframe_stopcodons)$lastEJ_dist)]
  
  # append 3' end of transcript to the first stop codon
  downUTRsize = length(thisqueryseq) - end(inframe_stopcodons[1])
  setORF = resizeTranscripts(queryCDS, end = downUTRsize) 
  
  ifelse(full.output == TRUE, 
         return(list(stopcodons = inframe_stopcodons, ORF = setORF)),
         return(list(stopcodons = inframe_stopcodons))
         )
}


#' Resize GRanges transcript object
#' 
#' @description 
#' This function will append/extend the 5' and 3' ends of a GRanges object which describes the 
#' exon ranges of a transcript or a CDS
#'
#' @usage resizeTranscripts(x, start = 0, end = 0)
#'
#' @param x GRanges object containing exon coordinates of a transcript or CDS
#' @param start Length of 5' end to truncate (positive val) or extend (negative val)
#' @param end Length of 3' end to truncate (positive val) or extend (negative val)
#'
#' @return a new GRanges transcript object 
#' @export
#'
#' @examples
#' resizeTranscripts(ptbp2_testData$exons$ENSMUST00000197833, start = 100)
#' 
resizeTranscripts <- function(x, start = 0, end = 0) {
  
  # get transcript strand information
  strand = as.character(strand(x))[1]
  exonsize = width(x)
  firstlastexons = c(1, length(x))
  mid = floor(length(x)/2)
  
  loop = TRUE
  while (loop == TRUE) {
    exonsize[firstlastexons] = exonsize[firstlastexons] - c(start, end)
    if (all(exonsize[firstlastexons] > 0)) {
      loop = FALSE
    } else {
      newval = ifelse(exonsize[firstlastexons] < 0, abs(exonsize[firstlastexons]), 0)
      start = newval[1]
      end = newval[2]
      exonsize = ifelse(exonsize > 0, exonsize, 0)
      firstlastexons[1] = ifelse(exonsize[firstlastexons[1]] > 0, firstlastexons[1], firstlastexons[1] + 1)
      firstlastexons[2] = ifelse(exonsize[firstlastexons[2]] > 0, firstlastexons[2], firstlastexons[2] - 1)
    }
  }
  
  #resize GRanges
  x[1:mid] = resize(x[1:mid], 
                    exonsize[1:mid], 
                    fix = ifelse(strand == '-', "start", "end"), 
                    ignore.strand = TRUE)
  x[mid:length(x)] = resize(x[mid:length(x)], 
                            exonsize[mid:length(x)], 
                            fix = ifelse(strand == '-', "end", "start"), 
                            ignore.strand = TRUE)
  
  x = x[width(x) > 0]
  return(x)
}