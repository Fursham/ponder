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
  diffSegments = suppressWarnings(indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")]))
  
  # test if the transcript contain same start codon as in CDS
    
    # do not test transcript and CDS do not overlap at all
  if (is.null(diffSegments)) {
    output = modifyList(output, 
                        list(txrevise_out = NA, 
                             annotatedStart = FALSE))
  } else if (length(diffSegments$refTx$upstream[diffSegments$refTx$upstream == TRUE]) == 0) {
    # transcripts that contain same start codon as CDS will return true for the above if statement
      
      # attempt to resize 5' end of transcritpt to the start codon
      # and also calculate size of first coding exon
    lenfirstsharedexon = width(diffSegments$shared_exons[1])
    upUTRsize = sum(width(diffSegments$testTx[diffSegments$testTx$upstream == TRUE]))
    setORF = resizeTranscripts(queryTX, start = upUTRsize)
    output = modifyList(output, 
                        list(txrevise_out = diffSegments, 
                             ORF = setORF, 
                             annotatedStart = TRUE, 
                             firstexlen = lenfirstsharedexon))
    } else {
      # transcripts that overlap with CDS but do not contain an annotated start codon 
      output = modifyList(output, 
                          list(txrevise_out = diffSegments, 
                               annotatedStart = FALSE))
  }
  ifelse(full.output == TRUE, return(output), return(output["annotatedStart"]))
}


#' Reconstruct CDS with alternate 5' exons
#'
#' @description This function will reconstruct a CDS sequence with a 5' end 
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
  
  reconstructedTx = sort(reduce(unlist(append(
    diffSegments$shared_exons, c(
    reduce(diffSegments[[1]][diffSegments[[1]]$upstream != TRUE]),
    reduce(diffSegments[[2]][diffSegments[[2]]$upstream == TRUE]))))),
    decreasing = queryStrand == '-')
  
  # prepares sequence of query
  thisqueryseq = unlist(Biostrings::getSeq(refsequence, reconstructedTx))
  
  # find all in-frame start codons in reconstructed transcript
    # it is vital that the last codon is a stop codon
  StartCodons = Biostrings::matchPattern(Biostrings::DNAString("ATG"), thisqueryseq)
  inFrameStartCodons = StartCodons[(length(reconstructedTx) - (end(StartCodons))) %%3 == 0,]
  
  # return ORF only if an in-frame start codon is found
  if (length(inFrameStartCodons) > 0) {
    
    # append 5' end of reconstructed transcript to the start codon
    upUTRsize = start(inFrameStartCodons[1]) - 1
    setORF = resizeTranscripts(reconstructedTx, start = upUTRsize)
    
    # obtain txrevise output
    combinedList = list(refTx = setORF, testTx = queryTx)
    diffSegments = indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
    out = list(ORF = setORF, txrevise_out = diffSegments, predictedStart = TRUE)
    
    if (is.null(diffSegments)) {
      out = list(ORF = NA, txrevise_out = NA, predictedStart = FALSE)
    }
    
  } else {
    out = list(ORF = NA, txrevise_out = NA, predictedStart = FALSE)
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
  
  # return NA if there is no unique internal and downstream alternative segments
  # else, construct new CDS with insertion of segments from query transcript

  if (length(diffSegments[[1]]$contained[diffSegments[[1]]$contained == TRUE]) == 0 &
      length(diffSegments[[1]]$downstream[diffSegments[[1]]$downstream == TRUE]) == 0 &
      length(diffSegments[[2]]$contained[diffSegments[[2]]$contained == TRUE]) == 0
      ) {
    Alternative_tx = FALSE
    augmentedCDS = sort(reduce(
      diffSegments$shared_exons),
      decreasing = as.character(strand(diffSegments$shared_exons))[1] == '-')
    
    ## WARNING
  } else {
    Alternative_tx = TRUE
    augmentedCDS = sort(reduce(unlist(append(
      diffSegments$shared_exons, 
      reduce(diffSegments[[2]][diffSegments[[2]]$upstream != TRUE])))),
      decreasing = as.character(strand(diffSegments$shared_exons))[1] == '-')
  }
    
    
    return(list(ORF = augmentedCDS, Alt_tx = Alternative_tx))
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
testClassicalNMD <- function(queryCDS, fasta, minmPTCdist = 50, full.output = FALSE){
  
  # obtain critical information
  queryStrand = as.character(strand(queryCDS))[1]
  exon_boundaries = cumsum(width(queryCDS))
  thisqueryseq = unlist(Biostrings::getSeq(fasta, queryCDS))
  
  # prepare a dict of stop codons for pattern matching
  list_stopcodons = Biostrings::DNAStringSet(c("TAA", "TAG", "TGA"))
  pdict_stopcodons = Biostrings::PDict(list_stopcodons)
  
  # search for in-frame stop codons
  allmatches = Biostrings::matchPDict(pdict_stopcodons, thisqueryseq) # check whether this list is named
  combinedmatches = unlist(allmatches)
  inframe_stopcodons = combinedmatches[end(combinedmatches) %% 3 == 0,]
  
  # calculate distance of stop codon to last exon junction
    # will return 0 if transcript is a single exon CDS
  lastEJ = ifelse(length(exon_boundaries) > 1,head(tail(exon_boundaries, n=2), n=1), 0)
  dist_stop_to_lastEJ = end(inframe_stopcodons) - lastEJ
  
  # update GRanges object 'inframe_stopcodons' with more information
  metadata = dplyr::data_frame(lastEJ_dist = dist_stop_to_lastEJ) %>%
    dplyr::mutate(is_NMD = ifelse(lastEJ_dist < -minmPTCdist, TRUE, FALSE)) %>%
    as.data.frame()
  elementMetadata(inframe_stopcodons) = metadata
  inframe_stopcodons = inframe_stopcodons[order(elementMetadata(inframe_stopcodons)$lastEJ_dist)]
  
  # append 3' end of transcript to the first stop codon
  downUTRsize = ifelse(length(inframe_stopcodons) > 0, length(thisqueryseq) - end(inframe_stopcodons[1]), 0)
  setORF = resizeTranscripts(queryCDS, end = downUTRsize) 
  
  # return output
  ifelse(full.output == TRUE, 
         return(list(stopcodons = inframe_stopcodons, ORF = setORF)),
         return(list(stopcodons = inframe_stopcodons))
         )
}





#' Title
#'
#' @param testCDS 
#' @param testTranscript 
#' @param fasta 
#'
#' @return
#' @export
#'
#' @examples
testNonClassicalNMD <- function(testCDS, testTranscript, fasta) {
  
  # test for non-coding segments
  combinedList = list(CDS = testCDS, Tx = testTranscript)
  nonCodingSegments = indentifyAddedRemovedRegions("CDS", "Tx", combinedList[c("CDS", "Tx")])
  
  # obtain size of 3'UTR
  threeUTR = sum(width(nonCodingSegments$Tx[nonCodingSegments$Tx$downstream == TRUE]))
  
  # test for uORF on 5'UTR
    # get sequence of 5'UTR
  strand = as.character(strand(testCDS))[1]
  fiveUTRGRanges = sort(nonCodingSegments$Tx[nonCodingSegments$Tx$upstream == TRUE], 
                        decreasing = strand == '-')
  fiveUTRseq = unlist(Biostrings::getSeq(fasta, fiveUTRGRanges))
  
    # prepare a dict of start and stop codons for pattern matching
  list_startstopcodons = Biostrings::DNAStringSet(c("ATG","TAA", "TAG", "TGA"))
  pdict_startstopcodons = Biostrings::PDict(list_startstopcodons)
  
    # search for start and stop codons
  allmatches = Biostrings::matchPDict(pdict_startstopcodons, fiveUTRseq)
  
    # this part attempts to add metadata to allmatches IRanges
  type = c(rep('Start', length(allmatches[[1]])), rep('Stop', length(unlist(allmatches))- length(allmatches[[1]])))
  combinedmatches = unlist(allmatches)
  elementMetadata(combinedmatches)$type = type
  frame = dplyr::data_frame(frame = end(combinedmatches)%%3) %>% as.data.frame()
  elementMetadata(combinedmatches)$frame = frame
  
  # return false if no start or stop codons are found
  if (length(combinedmatches) == 0) {
    return(list(uORF = FALSE, threeUTR = threeUTR))
  } else {
    # or else, test if subsets of in-frame Start and Stop codons are adjacent to each other
    uORF = NA
    for (i in 0:2) {
      ORForder = sort(combinedmatches[elementMetadata(combinedmatches)$frame == i])
      if (grepl('StartStop', paste(elementMetadata(ORForder)$type, collapse = ''))) {
        
        startGRanges = ORForder[elementMetadata(ORForder)$type == 'Start'][1]
        startIndex = match(startGRanges, ORForder)
        start = start(startGRanges) - 1
        
        newORF = ORForder[startIndex:length(ORForder)]
        end = length(fiveUTRseq) - end(newORF[elementMetadata(newORF)$type == 'Stop'][1])
        
        ORFcoord = resizeTranscripts(fiveUTRGRanges, start, end)
        
        uORF = paste(ranges(ORFcoord), collapse = ';')
        break
      }
    }
  }
  
  return(list(uORF = uORF, threeUTR = threeUTR))
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
  startend = c(start,end)
  
  for (i in 1:2) {
    loop = TRUE
    while (loop == TRUE) {
      exonsize[firstlastexons[[i]]] = exonsize[firstlastexons[[i]]] - startend[[i]]
      if (all(exonsize >= 0)) {
        loop = FALSE
      } else {
        newval = ifelse(exonsize[firstlastexons[[i]]] < 0, abs(exonsize[firstlastexons[[i]]]), 0)
        startend[[i]] = newval
        exonsize = ifelse(exonsize > 0, exonsize, 0)
        if (i == 1) {
          firstlastexons[1] = ifelse(exonsize[firstlastexons[1]] > 0, firstlastexons[1], firstlastexons[1] + 1)
        } else if (i == 2) {
          firstlastexons[2] = ifelse(exonsize[firstlastexons[2]] > 0, firstlastexons[2], firstlastexons[2] - 1)
        }
      }
    
    }
    if (i == 1) {
      x[1:length(x)] = resize(x[1:length(x)], 
                        exonsize[1:length(x)], 
                        fix = ifelse(strand == '-', "start", "end"), 
                        ignore.strand = TRUE)
    } else if (i == 2) {
      x[1:length(x)] = resize(x[1:length(x)], 
                                exonsize[1:length(x)], 
                                fix = ifelse(strand == '-', "end", "start"), 
                                ignore.strand = TRUE)
    }
  }
  # return exons with width of more than 0
  x = x[width(x) > 0]
  return(x)
}





classifyAltSegments <- function(transcript1, transcript2, txrevise_out) {
  
  # check if txrevise_out input is provided, and build one if not.
  if (missing(txrevise_out)) {
    if (missing(transcript1) | missing(transcript2)) {
      stop('Please provide input GRanges objects')
    }
    combinedList = list(transcript1 = transcript1, transcript2 = transcript2)
    diffSegments = indentifyAddedRemovedRegions("transcript1", "transcript2", combinedList[c("transcript1", "transcript2")])
  } else {
    diffSegments = txrevise_out
  }
  
  # escape if both transcripts do not have overlapping exons
  if (is.null(diffSegments)) {
    return(NULL) # warning too maybe?
  }
  
  # get coverage of overlaps. size of shared exons / size of combined exons
  combinedLength = sum(width(reduce(append(transcript1, transcript2))))
  sharedLength = sum(width(diffSegments[[3]]))
  coverage = round((sharedLength/combinedLength), 3)
  diffSegments = c(diffSegments, Coverage = coverage)
  
  
  # cycle transcript specific segments
  for (i in 1:2) {
    AS_class = c()
    
    # skip transcript with no alternative segments
    if (length(diffSegments[[i]]) == 0) {
      next
    } else {
      for (j in 1:length(diffSegments[[i]])) {
        thisexon = diffSegments[[i]][j]
        
        # test for alternative first and last exons first, then internal exons
        if (elementMetadata(thisexon)$upstream == TRUE) {
          AS_class = c(AS_class, 'AF')
          next
        } else if (elementMetadata(thisexon)$downstream == TRUE) {
          AS_class = c(AS_class, 'AL')
        } else {
          
          # test for alt splice sites
          mergedsegment = reduce(append(ranges(thisexon), ranges(diffSegments[[3]])))
          if (length(mergedsegment) == length(diffSegments[[3]])) {
            # alternative splice sites
            
            # test alternative 5'
            if (any(!countOverlaps(mergedsegment, ranges(diffSegments[[3]]), type = 'start'))) {
              AS_class = c(AS_class, 'A5')
              next
            } 
            # test alternative 3'
            else if (any(!countOverlaps(mergedsegment, ranges(diffSegments[[3]]), type = 'end'))) {
              AS_class = c(AS_class, 'A3')
              next
            }
            
          } else if (length(mergedsegment) > length(diffSegments[[3]])){
            # test casette vs MX
            
            # get flanking exons
            sortedmergedsegment = sort(mergedsegment, decreasing = as.character(strand(thisexon)) == '-')
            exonindex = match(ranges(thisexon), sortedmergedsegment)
            flankexonsindex = c(exonindex-1, exonindex+1)
            
            # 
            flankcoords = c(sortedmergedsegment[flankexonsindex[1]], sortedmergedsegment[flankexonsindex[2]])
            mergedflanks = range(flankcoords)
            
            if (i == 1) {
              othertx = diffSegments[[2]]
            } else {
              othertx = diffSegments[[1]]
            }

            if (length(othertx) == 0) {
              AS_class = c(AS_class, 'CE')
              next
            }
            
            matchedindex = countOverlaps(mergedflanks, ranges(othertx))
            matchedsegments = othertx[matchedindex]
            
            # if there are other segments in the other transcript, this could be an MX
            if (length(matchedsegments) > 0) {
              tag = FALSE
              for (k in 1:length(matchedsegments)) {
                # test for possibility of the other segment to be a5' or a3'
                newmergedsegment = reduce(append(ranges(matchedsegments[k]), flankcoords))
                
                # if newmergedsegment > 2, it means that the segment is not a5' or a3'
                if (length(newmergedsegment) > 2) {
                  AS_class = c(AS_class, 'MX')
                  tag = TRUE
                  break
                }
              }
              # if none of the matchedsegments are MX, tag exon as 'CE'
              if (tag == FALSE) {
                AS_class = c(AS_class, 'CE')
                next
              }
            }
            # if there are no other segments in other transcript, it's most prob a casette exon
            else {
              AS_class = c(AS_class, 'CE')
              next
            }
          }
          # if length of mergedsegment is smaller, it's most likley and intron retention
          else if (length(mergedsegment) < length(diffSegments[[3]])){
            AS_class = c(AS_class, 'IR')
            next
          }
        }
      }
    }
    
  if (i == 1) {
    AS_class = tolower(AS_class)
  } else {
    AS_class = toupper(AS_class)
  }
  elementMetadata(diffSegments[[i]]) = cbind(AS_class)
  }
  return(diffSegments)
}



