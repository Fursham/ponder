#' Reconstruct CDS with alternative internal and downstream segments
#'
#' @description 
#' This function will add alternative segments from a query transcript into 
#' a reference CDS from the same gene and generate a new ORF
#' 
#' Note: This function will not insert/remove first exons. 
#' 
#' @param queryTranscript A GRanges object containing exon structure of a query transcript from
#' the same gene family
#' @param refCDS A GRanges object containing CCDS information of a gene family
#' @param fasta Fasta sequence of the genome
#' @param txrevise_out 
#' Optional. 
#' GRangesobject from 'indentifyAddedRemovedRegions()' 
#' or 'testTXforStart()' output. ?indentifyAddedRemovedRegions or ?testTXforStart
#' for more information.
#' Arguments refCDS and queryTranscript is not mandatory if this argument is provided
#' 
#' 
#' @return
#' A list containing:
#' (1) A GRanges Object of new ORF, or NA if no ORF is found
#' (2) TRUE/FALSE object on whether ORF is an alternative CDS transcript

#' @export
#' @import Biostrings
#' 
#' @author Fursham Hamid
#'
#' @examples
#' 
#' library("BSgenome.Mmusculus.UCSC.mm10")
#' reconstructCDS(ptbp2Data$transcripts$ENSMUST00000197833, ptbp2Data$refCDS, fasta = BSgenome.Mmusculus.UCSC.mm10)
#' 
#' 
reconstructCDS <- function(queryTranscript, refCDS, fasta, txrevise_out = NULL, gene_id, transcript_id){
  
  # check if txrevise_out input is provided, and build one if not.
  if (is.null(txrevise_out)) {
    if (missing(refCDS) | missing(queryTranscript)) {
      stop('Please provide input GRanges objects')
    }
    combinedList = list(refTx = refCDS, testTx = queryTranscript)
    diffSegments = indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
  } else {
    diffSegments = txrevise_out
  }
  
  # prepare output
  output = list(ORF_considered = NA, Alt_tx = NA)
  
  # return NA if there is no unique internal and downstream alternative segments
  # else, construct new CDS with insertion of segments from query transcript
  
  # Basically, this part tests if there is any unique internal and downstream segments
  if (length(diffSegments[[1]]$contained[diffSegments[[1]]$contained == TRUE]) == 0 &
      length(diffSegments[[1]]$downstream[diffSegments[[1]]$downstream == TRUE]) == 0 &
      length(diffSegments[[2]]$contained[diffSegments[[2]]$contained == TRUE]) == 0) {
    
    # if not, it will return the refCDS as the CDS and return alternative_tx as false
    Alternative_tx = FALSE
    augmentedCDS = sort(reduce(
      diffSegments$shared_exons),
      decreasing = as.character(strand(diffSegments$shared_exons))[1] == '-')
    
    ## POSSIBLE WARNING
  } else {
    # if there are, construct a new exon structure with removal/addition of the alternative segments 
    Alternative_tx = TRUE
    augmentedCDS = sort(append(
      diffSegments$shared_exons, 
      reduce(diffSegments[[2]][diffSegments[[2]]$upstream != TRUE])),
      decreasing = as.character(strand(diffSegments$shared_exons))[1] == '-')
    
    # this part will correct the open reading frame
      # obtain critical information on the new CDS
    queryStrand = as.character(strand(augmentedCDS))[1]
    thisqueryseq = unlist(Biostrings::getSeq(fasta, augmentedCDS))
    
      # prepare a dict of stop codons for pattern matching
    list_stopcodons = Biostrings::DNAStringSet(c("TAA", "TAG", "TGA"))
    pdict_stopcodons = Biostrings::PDict(list_stopcodons)
    
      # search for in-frame stop codons
    allmatches = Biostrings::matchPDict(pdict_stopcodons, thisqueryseq)
    combinedmatches = unlist(allmatches)
    inframe_stopcodons = sort(combinedmatches[end(combinedmatches) %% 3 == 0,])
    
      # append 3' end of transcript to the first stop codon if found
    if (length(inframe_stopcodons) > 0) {
      downUTRsize = length(thisqueryseq) - end(inframe_stopcodons[1])
      augmentedCDS = resizeTranscripts(augmentedCDS, end = downUTRsize)
      augmentedCDS = augmentedCDS %>% as.data.frame() %>%
        dplyr::mutate(type = 'CDS', gene_id = gene_id, transcript_id = transcript_id) %>%
        dplyr::mutate(phase = cumsum(width%%3)%%3)
      augmentedCDS$phase = c(0, head(augmentedCDS$phase, - 1))
      augmentedCDS = makeGRangesFromDataFrame(augmentedCDS, keep.extra.columns = TRUE)

    } else {
      # if a stop codon is not found, return NA
      augmentedCDS = NA
    }
  }
  output = modifyList(output, 
                      list(ORF_considered = augmentedCDS, Alt_tx = Alternative_tx))
  return(output)
}





#' Test ORF for NMD
#' 
#' @description 
#' This function will test for NMD features of a transcript with a given open reading frame
#' 
#'
#' @param queryCDS GRanges Object containing CDS structure of a transcript
#' @param queryTranscript GRanges Object containing transcript structure of the same transcript
#' @param distance_stop_EJ NMD Threshold. Distance of stop codon to last exon-exon junction. Default: 50
#' @param other_features TRUE/FALSE argument to report other NMD features. Default: FALSE
#' @param fasta Fasta sequence file. Mandatory if other_features == TRUE
#' 
#' @return 
#' A list containing information on: 
#' (1) NMD nature of transcript
#' (2) Distance of stop codon to last exon-exon junction
#' if other features of NMD is to be returned:
#' (3) coordinates of uORF
#' (4) coordinates of uATG
#' (5) frame of uATG
#' 
#' @import Biostrings
#' @author Fursham Hamid
#' @export
#' 
#' @examples
#' 
#' library("BSgenome.Mmusculus.UCSC.mm10")
#' testNMD(ptbp2Data$refCDS, ptbp2Data$transcripts$ENSMUST00000029780)
#' testNMD(ptbp2Data$refCDS, ptbp2Data$transcripts$ENSMUST00000029780, other_features = TRUE, fasta = BSgenome.Mmusculus.UCSC.mm10)
#' 
#' testNMD(ptbp2Data$skipE10CDS, ptbp2Data$transcripts$ENSMUST00000197833)
#' testNMD(ptbp2Data$skipE10CDS, ptbp2Data$transcripts$ENSMUST00000197833, other_features = TRUE, fasta = BSgenome.Mmusculus.UCSC.mm10)
#' 
#' 
testNMD <- function(queryCDS, queryTranscript, distance_stop_EJ = 50, other_features = FALSE, fasta){
  
  # prepare output list
  output = list(is_NMD = as.logical(NA), dist_to_lastEJ = as.numeric(NA))
  if (other_features == TRUE) {
    output = modifyList(output, list(uORF = as.character(NA), 
                                     threeUTR = as.numeric(NA), 
                                     uATG = as.character(NA), 
                                     uATG_frame = as.character(NA)))
  }
  
  if (is.na(queryCDS[1])) {
    return(output)
  }
  
  # sort queryCDS, by exon order (just in case)
  strand = as.character(strand(queryCDS))[1]
  queryCDS = sort(queryCDS, decreasing = strand == '-')
  queryTranscript = sort(queryTranscript, decreasing = strand == '-')
  
  # get noncoding segments from indentifyAddedRemovedRegions
  combinedList = list(CDS = queryCDS, Tx = queryTranscript)
  noncodingSegments = indentifyAddedRemovedRegions("CDS", "Tx", combinedList[c("CDS", "Tx")])
  
  # obtain coordinates of CDS and its 3'UTR
  cds3UTR = sort(append(
    noncodingSegments$shared_exons, 
    reduce(noncodingSegments$Tx[noncodingSegments$Tx$downstream == TRUE])),
    decreasing = (strand == '-'))
  
  # obtain distance of stop codon to last exon-exon junction and modify output list
  stop_codon_position = sum(width(queryCDS))
  exon_boundaries = cumsum(width(cds3UTR))
  lastEJ = ifelse(length(exon_boundaries) > 1,head(tail(exon_boundaries, n=2), n=1), 0)
  dist_stop_to_lastEJ = stop_codon_position - lastEJ
  is_NMD = ifelse(dist_stop_to_lastEJ < -(distance_stop_EJ), TRUE, FALSE)
  output = modifyList(output, list(is_NMD = is_NMD, dist_to_lastEJ = dist_stop_to_lastEJ))
  
  # test for other features if activated
  if (other_features == TRUE) {
    if (missing(fasta)) {
      stop('Please provide fasta sequence')
    }
    # obtain size of 3'UTR
    threeUTR = sum(width(noncodingSegments$Tx[noncodingSegments$Tx$downstream == TRUE]))
    
    # test for uORF on 5'UTR
      # get sequence of 5'UTR
    fiveUTRGRanges = sort(noncodingSegments$Tx[noncodingSegments$Tx$upstream == TRUE], 
                          decreasing = strand == '-')
    list_startstopcodons = Biostrings::DNAStringSet(c("ATG","TAA", "TAG", "TGA"))
    pdict_startstopcodons = Biostrings::PDict(list_startstopcodons)
    
    # this part will test the presence of uORFs and uATGs in the 5UTR
    
      # prepare output variables
    uORF = as.character(NA)
    uATG = as.character(NA)
    uATG_frame = as.character(NA)
    
    # cycle each frame to find 5UTR and uATGs
    for (i in 0:2) {
      thisfiveUTRGRanges = fiveUTRGRanges
      # get sequence of 5UTR and search for all start and stop codons
      fiveUTRseq = unlist(Biostrings::getSeq(fasta, thisfiveUTRGRanges))
      allmatches = Biostrings::matchPDict(pdict_startstopcodons, fiveUTRseq)
      
      # this part adds type and frame information into allmatches as a metadata
      type = c(rep('Start', length(allmatches[[1]])), rep('Stop', length(unlist(allmatches))- length(allmatches[[1]])))
      combinedmatches = unlist(allmatches)
      elementMetadata(combinedmatches)$type = type
      frame = dplyr::data_frame(frame = (length(fiveUTRseq) - end(combinedmatches))%%3) %>% as.data.frame()
      elementMetadata(combinedmatches)$frame = frame
      
      # break loop if there are no start and stop codons for frame i
      if ((length(combinedmatches) == 0)){
        next
      } else if(!any(mcols(combinedmatches)$frame == i)) {
        next
      }
      
      shiftype = c('Stop', head(type, length(type)-1))
      elementMetadata(combinedmatches)$shiftype = shiftype
      ORForder = sort(combinedmatches[elementMetadata(combinedmatches)$frame == i])
      if ((length(ORForder[elementMetadata(ORForder)$type == 'Start']) == 0)) {
        next
      } 
      ORForder = ORForder[elementMetadata(ORForder)$type != elementMetadata(ORForder)$shiftype]
      if (length(ORForder) == 0) {
        next
      }
      
      for (j in seq(1, length(ORForder), 2)) {
        
        startGRanges = ORForder[j]
        startFrame = elementMetadata(startGRanges)$frame
        startlength = start(startGRanges) - 1
        
        if (match(startGRanges, ORForder) == length(ORForder)) {
          # return uATG
          endlength = length(fiveUTRseq) - end(startGRanges)
          startCoords = resizeTranscripts(thisfiveUTRGRanges, startlength, endlength)
          uATG = c(uATG, paste(ranges(startCoords)))
          uATG_frame = c(uATG_frame, startFrame)
        } else {
          # return uORF
          stopGRanges = ORForder[(j+1)]
          endlength = length(fiveUTRseq) - end(stopGRanges)
          ORFcoord = resizeTranscripts(thisfiveUTRGRanges, startlength, endlength)
          uORF = c(uORF, paste(ranges(ORFcoord), collapse = ';')) # need to remove NA first later
        }
      }
    }
    # update output list
    uORF = ifelse(all(is.na(uORF)), FALSE, paste(uORF[!is.na(uORF)], collapse = '|'))
    uATG = ifelse(all(is.na(uATG)), FALSE, paste(uATG[!is.na(uATG)], collapse = '|'))
    uATG_frame = ifelse(all(is.na(uATG_frame)), FALSE, paste(uATG_frame[!is.na(uATG_frame)], collapse = '|'))
    
    output = modifyList(output, list(uORF = as.character(uORF), 
                                     threeUTR = as.numeric(threeUTR), 
                                     uATG = as.character(uATG), 
                                     uATG_frame = as.character(uATG_frame)))
  }
  # return output
  return(output)
}
    
      


#' Classify alternative splicing features
#' 
#' @description 
#' This function is an extension of the function indentifyAddedRemovedRegions by classifying
#' of the alternative segments
#'
#' @param transcript1 GRanges Object containing exon structure of transcript 1
#' @param transcript2 GRanges Object containing exon structure of transcript 2
#' @param txrevise_out Optional.
#'
#' @return a list of GRanges objects similar to 'indentifyAddedRemovedRegions()' output (txrevise_out). 
#' ?indentifyAddedRemovedRegions for more information
#' @export
#'
#' @examples
#' 
#' classifyAltSegments(ptbp2Data$transcripts$ENSMUST00000029780, ptbp2Data$transcripts$ENSMUST00000197833)
#' 
classifyAltSegments <- function(transcript1, transcript2, txrevise_out = NULL) {
  
  # check if txrevise_out input is provided, and build one if not.
  if (is.null(txrevise_out)) {
    if (missing(transcript1) | missing(transcript2)) {
      stop('Please provide input GRanges objects')
    }
    combinedList = list(transcript1 = transcript1, transcript2 = transcript2)
    diffSegments = indentifyAddedRemovedRegions("transcript1", "transcript2", combinedList[c("transcript1", "transcript2")])
  } else {
    diffSegments = txrevise_out
  }
  
  # return if both transcripts do not have overlapping exons
  if (is.null(diffSegments)) {
    return(NULL) # warning too maybe?
  }
  
  # get coverage length of overlaps
  coverage = sum(width(diffSegments[[3]])) / sum(width(c(ranges(diffSegments[[1]]), ranges(diffSegments[[3]]))))
  diffSegments = c(diffSegments, Shared_coverage = coverage)
  
  # test for alternate segments in each transcripts, one at a time
  strand = as.character(strand(transcript1)[1])
  totest = sapply(list(diffSegments[[1]], diffSegments[[2]]), function(x) {
    ifelse(length(x) == 0, FALSE, TRUE)
  })
  for (i in (1:2)[totest]) {
    
    # set a vector to store alternative splicing classes
    AS_class = apply(as.data.frame(diffSegments[[i]]), 1, function(x) {
      
      # convert row to list for easy referencing
      thisexon = x %>% as.list()
      thisGRanges = GRanges(seqnames = thisexon$seqnames, 
                            ranges = IRanges(as.numeric(thisexon$start), as.numeric(thisexon$end)),
                            strand = thisexon$strand)
      
      # test for alternative first and last exons first, then internal exons
      if (as.numeric(thisexon$upstream) == 1) {
        
        # alternative segment could be an alternate transcription start or alternate first exon
        #   to distinguish this, we merge the exon with the shared exon and determine number of exons combined
        mergedsegment = reduce(append(ranges(thisGRanges), ranges(diffSegments[[3]])))
        if (length(mergedsegment) == length(diffSegments[[3]])) {
          exonindex = which(countOverlaps(combinedList[[i]], thisGRanges) == 1)
          if (exonindex == 1) {
              return('ATS')
            } else {
              return('AF')
            }
        } else {
          return('AF')
        }
      } else if (as.numeric(thisexon$downstream) == 1) {
        mergedsegment = reduce(append(ranges(thisGRanges), ranges(diffSegments[[3]])))
        if (length(mergedsegment) == length(diffSegments[[3]])) {
          exonindex = which(countOverlaps(combinedList[[i]], thisGRanges) == 1)
          totalexonnum = length(combinedList[[i]])
          
          if (exonindex == totalexonnum) {
              return('APA')
            } else {
              return('AL')
            }
        } else {
          return('AL')
        }
      } else {
        
        # classify alternative internal segments
        mergedsegment = reduce(append(ranges(thisGRanges), ranges(diffSegments[[3]])))
        
        # idea here is that after merging the alternate segment with shared exons,
        #   if the number of exons remain the same, alternate segment is an alternative splice sites
        #   if the number of exons increase by 1, segment is either casette or mutually exclusive
        #   if the number of exons decrease by 1, segment is a retained intron
        
        # alternative splice sites
        if (length(mergedsegment) == length(diffSegments[[3]])) {
          
          # to distinguish alternative 5' and 3' splice site events, we determine if there is an 
          # overlap betweeen mergedsegment and shared_exons at either the start or end of the exons
          
          if (any(!countOverlaps(mergedsegment, ranges(diffSegments[[3]]), type = 'start'))) {
            if (strand == '-') {
              return('A5')
            } else {
              return('A3')
            }
          } 
          else if (any(!countOverlaps(mergedsegment, ranges(diffSegments[[3]]), type = 'end'))) {
            if (strand == '-') {
              return('A3')
            } else {
              return('A5')
            }
          }
          
          # casette or mutually exclusive events 
        } else if (length(mergedsegment) > length(diffSegments[[3]])){
          
          # get indexes and coordinates of flanking exons
          sortedmergedsegment = sort(mergedsegment, decreasing = as.character(strand(thisGRanges)) == '-')
          exonindex = match(ranges(thisGRanges), sortedmergedsegment)
          flankexonsindex = c(exonindex-1, exonindex+1)
          flankcoords = c(sortedmergedsegment[flankexonsindex[1]], sortedmergedsegment[flankexonsindex[2]])
          mergedflanks = range(flankcoords)
          
          # get the list of unique segments from other transcript
          if (i == 1) {
            othertx = diffSegments[[2]]
          } else {
            othertx = diffSegments[[1]]
          }
          
          # if there no other unique segments, return as casette exon
          if (length(othertx) == 0) {
            return('CE')
          }
          
          # find all the unique segments that fall within the flanking exons
          matchedindex = countOverlaps(mergedflanks, ranges(othertx))
          matchedsegments = othertx[matchedindex]
          
          # if there are other segments in the other transcript, this could be an MX
          if (length(matchedsegments) > 0) {

            testMX = apply(as.data.frame(matchedsegments), 1, function(y) {
              # test for possibility of the other segment to be a5' or a3'
              y = as.list(y)
              newmergedsegment = reduce(append(IRanges(as.numeric(y$start), as.numeric(y$end)), flankcoords))
              
              # if newmergedsegment > 2, it means that the segment is not a5' or a3'
              if (length(newmergedsegment) > 2) {
                return('MX')
              } else {
                return(NA)
              }
            })
            if ('MX'%in%testMX) {
              return('MX')
            } else {
              return('CE')
            }
          }
          # if there are no other segments in the flanking exons in other transcript, it's most prob a casette exon
          else {
            return('CE')
          }
        }
        # if length of mergedsegment is smaller, it's most likley and intron retention
        else if (length(mergedsegment) < length(diffSegments[[3]])){
          return('IR')
        }
      }
    })

    # Set classes in transcript 1 as lower case and transcript 2 as upprecase
    if (i == 1) {
      AS_class = tolower(AS_class)
    } else {
      AS_class = toupper(AS_class)
    }
    elementMetadata(diffSegments[[i]]) = cbind(AS_class)
  }
  return(diffSegments)
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
#' @import GenomicRanges
#' @return a new GRanges transcript object 
#' @export
#'
#' @examples
#' 
#' resizeTranscripts(ptbp2Data$transcripts$ENSMUST00000197833, start = 100)
#' 
resizeTranscripts <- function(x, start = 0, end = 0) {
  
  if (sum(width(x)) < (start + end)) {
    stop('Appending length is larger than size of transcript')
  }
  
  # retrieve important information
  strand = as.character(strand(x))[1]
  firstlastexons = c(1, length(x))  # initiate index of first and last exons
  startend = c(start,end) # store value for appending
  
  totest = sapply(startend, function(x) {
    ifelse(x > 0, TRUE, FALSE)
  })
  # append the start followed by the end
  for (i in (1:2)[totest]) {
    exonsize = width(x) # get the size of each exons as a list
    loop = TRUE
    while (loop == TRUE) {
      # change the size of the first or last exon
      exonsize[firstlastexons[[i]]] = exonsize[firstlastexons[[i]]] - startend[[i]]
      
      # if all of the exon sizes are more than 0, it means that the subtraction occurs within an exon
      #   if there are some exon with negative sizes, it means that the subtraction bleed into the next exon
      if (all(exonsize >= 0)) {
        loop = FALSE
      } else {
        
        # set the remaining start/end value to be subtracted
        newval = ifelse(exonsize[firstlastexons[[i]]] < 0, abs(exonsize[firstlastexons[[i]]]), 0)
        startend[[i]] = newval
        
        # convert negative widths to 0
        exonsize = ifelse(exonsize > 0, exonsize, 0)
        
        # and change index of exon to be subtracted
        if (i == 1) {
          firstlastexons[1] = ifelse(exonsize[firstlastexons[1]] > 0, firstlastexons[1], firstlastexons[1] + 1)
        } else if (i == 2) {
          firstlastexons[2] = ifelse(exonsize[firstlastexons[2]] > 0, firstlastexons[2], firstlastexons[2] - 1)
        }
      }
      
    }
    # resize the transcripts based on its new width
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





