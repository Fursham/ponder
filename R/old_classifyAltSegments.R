
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


