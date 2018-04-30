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
reconstructCDS <- function(queryTranscript, refCDS, fasta, txrevise_out = NULL){
  
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
    augmentedCDS = sort(reduce(unlist(append(
      diffSegments$shared_exons, 
      reduce(diffSegments[[2]][diffSegments[[2]]$upstream != TRUE])))),
      decreasing = as.character(strand(diffSegments$shared_exons))[1] == '-')
    
    # this part will correct the open reading frame
      # obtain critical information on the new CDS
    queryStrand = as.character(strand(augmentedCDS))[1]
    exon_boundaries = cumsum(width(augmentedCDS))
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
    } else {
      # if a stop codon is not found, return NA
      augmentedCDS = NA
    }
  }
  return(list(ORF_considered = augmentedCDS, Alt_tx = Alternative_tx))
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
  output = list()
  
  # sort queryCDS, by exon order (just in case)
  strand = as.character(strand(queryCDS))[1]
  queryCDS = sort(queryCDS, decreasing = strand == '-')
  queryTranscript = sort(queryTranscript, decreasing = strand == '-')
  
  # get noncoding segments from indentifyAddedRemovedRegions
  combinedList = list(CDS = queryCDS, Tx = queryTranscript)
  noncodingSegments = indentifyAddedRemovedRegions("CDS", "Tx", combinedList[c("CDS", "Tx")])
  
  # obtain coordinates of CDS and its 3'UTR
  cds3UTR = sort(reduce(unlist(append(
    noncodingSegments$shared_exons, 
    reduce(noncodingSegments$Tx[noncodingSegments$Tx$downstream == TRUE])))),
    decreasing = as.character(strand(noncodingSegments$shared_exons))[1] == '-')
  
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
    uORF = NA
    uATG = NA
    uATG_frame = NA
    
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
        startFrame = elementMetadata(startGRanges)$frame[[1]]
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
    
    output = modifyList(output, list(uORF = uORF, threeUTR = threeUTR, uATG = uATG, uATG_frame = uATG_frame))
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
  for (i in 1:2) {
    
    # set a vector to store alternative splicing classes
    AS_class = c()
    
    # skip transcript with no alternative segments
    if (length(diffSegments[[i]]) == 0) {
      next
    } else {
      for (j in 1:length(diffSegments[[i]])) {
        thisexon = diffSegments[[i]][j]
        
        # test for alternative first and last exons first, then internal exons
        if (elementMetadata(thisexon)$upstream == TRUE) {
          
          # alternative segment could be an alternate transcription start or alternate first exon
          #   to distinguish this, we merge the exon with the shared exon and determine number of exons combined
          mergedsegment = reduce(append(ranges(thisexon), ranges(diffSegments[[3]])))
          if (length(mergedsegment) == length(diffSegments[[3]])) {
            AS_class = c(AS_class, 'ATS')
            next
          } else {
            AS_class = c(AS_class, 'AF')
            next
          }
        } else if (elementMetadata(thisexon)$downstream == TRUE) {
          mergedsegment = reduce(append(ranges(thisexon), ranges(diffSegments[[3]])))
          if (length(mergedsegment) == length(diffSegments[[3]])) {
            AS_class = c(AS_class, 'APA')
            next
          } else {
            AS_class = c(AS_class, 'AL')
            next
          }
        } else {
          
          # classify alternative internal segments
          mergedsegment = reduce(append(ranges(thisexon), ranges(diffSegments[[3]])))
          
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
                AS_class = c(AS_class, 'A5')
              } else {
                AS_class = c(AS_class, 'A3')
              }
              next
            } 
            else if (any(!countOverlaps(mergedsegment, ranges(diffSegments[[3]]), type = 'end'))) {
              if (strand == '-') {
                AS_class = c(AS_class, 'A3')
              } else {
                AS_class = c(AS_class, 'A5')
              }
              next
            }
            
          # casette or mutually exclusive events 
          } else if (length(mergedsegment) > length(diffSegments[[3]])){
            
            # get indexes and coordinates of flanking exons
            sortedmergedsegment = sort(mergedsegment, decreasing = as.character(strand(thisexon)) == '-')
            exonindex = match(ranges(thisexon), sortedmergedsegment)
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
              AS_class = c(AS_class, 'CE')
              next
            }
            
            # find all the unique segments that fall within the flanking exons
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
            # if there are no other segments in the flanking exons in other transcript, it's most prob a casette exon
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
  
  # retrieve important information
  strand = as.character(strand(x))[1]
  firstlastexons = c(1, length(x))  # initiate index of first and last exons
  startend = c(start,end) # store value for appending
  
  # append the start followed by the end
  for (i in 1:2) {
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








#' Match Gene IDs from query GTF/GFF3 file
#' 
#' @description 
#' This function will match and correct Gene IDs from a query assembled transcript file, using a 
#' transcript annotation as reference. 
#' 
#' The default approach to this correction relies on finding overlaps between transcripts in query with
#' transcripts in reference. Using this method alone could result in false positive matches (19 percent false positives).
#' To improve this, users have an option to activate two additional layers of matching.
#' (1) Matching by ENSEMBL Gene_IDs. If both query and reference transcript annotations containg Ensembl-style
#' Gene IDs, this program will try to match both IDs in a less stringent manner. This correction can be invoked
#' by providing the 'primary_gene_id' argument
#' 
#' (2) Matching by secondary Gene_IDs. Depending on the transcript assembly program, GTF/GFF3 annotations 
#' may contain additional comments on the transcript information. This may include a distinct 
#' secondary Gene ID annotation that potentially matches with the reference. To invoke this correction,
#' provide 'primary_gene_id' and 'secondary_gene_id' arguments. To determine if your transcript assembly contain 
#' possible secondary Gene IDs, try importing query GTF file using rtracklayer package and check its metadata 
#' columns
#' 
#' 
#' @param query 
#' Assembled transcripts in query. Argument can be provided as name of GTF/GFF3 file in current working directory,
#' or path_to_file. Argument can also be provided as a GenomicRanges object containing transcript and/or exon ranges.
#' @param ref 
#' Annotated reference transcripts. This program contain gencode basic annotations from mm10 and hg38 assemblies  
#' that will be loaded using 'gencode_basics[['mm10']] or gencode_basics[['hg38']] as input to argument. Argument can also be provided as name of 
#' GTF/GFF3 file in current working directory, or path_to_file.
#' @param query_format 
#' Format of query file (GTF/GFF3). If this argument is not provided and input to query argument
#' is a filename, program will attempt to read file extension.
#' @param ref_format 
#' Format of reference file (GTF/GFF3). If this argument is not provided and input to query argument
#' is a filename, program will attempt to read file extension.
#' @param primary_gene_id 
#' Name of the primary gene id in query file. Input to this argument is typically 'gene_id'
#' @param secondary_gene_id 
#' Name of the secondary gene id in query file. Example of input to this arguement is 'ref_gene_id'
#' @param makefile 
#' TRUE/FALSE input. If TRUE, program will output a GTF file. Else, program will return a GenomicRanges
#' object with matched Gene IDs
#' 
#' @param outputfile 
#' Name of output file (Default: 'matched_geneIDs.gtf'). File will be saved in current working directory
#'
#' @return
#' GTF file with matched Gene IDs if makefile == TRUE. 
#' GRanges object if makefile == FALSE.
#' 
#' @export
#'
#' @examples
#' 
#' matchGeneIDs(testData, gencode_basics[['mm10']], primary_gene_id = 'gene_id', secondary_gene_id = 'ref_gene_id', clusters = 4)
#' 
#' 
#' 
matchGeneIDs <- function(query, ref, query_format = NULL, ref_format=NULL, primary_gene_id=NULL, secondary_gene_id=NULL, makefile = TRUE, outputfile = 'matched_geneIDs.gtf', clusters = 4) {
  
  # check the nature of query and ref. if it's a GRanges, proceed with analysis, if not, attempt to import
  if (is(query, 'GenomicRanges')) {
    # return error if input and reference do not contain gene_id metadata
    if (is.null(elementMetadata(query)$gene_id)) {
      if (makefile == FALSE) {
        stopLog('Gene_ID metadata error: Please ensure input assembled transcripts contain gene_id metadata\n', logf)
      } else {
        stop('Gene_ID metadata error: Please ensure input assembled transcripts contain gene_id metadata\n')
      }
    }
    inputGRanges = query
  } else {
    # try to import query
    if (!file.exists(query)){
      if (makefile == FALSE) {
        stopLog('Input transcript file do not exist', logf)
      } else {
        stop('Input transcript file do not exist')
      }
    }
    
    # use file extension format if provided by user, else try to extract file format from filename
    if (!is.null(query_format)) {
      fileformat = query_format
    } else {
      fileformat = tail(unlist(strsplit(query, '\\.')), "1")
      if (!fileformat%in%c('gff3','gff','gtf','bed')) {
        if (makefile == FALSE) {
          stopLog('Incorrect input file extension format', logf)
        } else {
          stop('Incorrect input file extension format')
        }
      } 
    }
    inputGRanges = rtracklayer::import(query, format = fileformat)
  }
  
  if (is(ref, 'GenomicRanges')) {
    # return error if input and reference do not contain gene_id metadata
    if (is.null(elementMetadata(ref)$gene_id)) {
      if (makefile == FALSE) {
        stopLog('Gene_ID metadata error: Please ensure reference assembly contain gene_id metadata', logf)
      } else {
        stop('Gene_ID metadata error: Please ensure reference assembly contain gene_id metadata')
      }
    }
    basicGRanges = ref
  } else {
    # try to import ref
    if (!file.exists(ref)){
      if (makefile == FALSE) {
        stopLog('Reference transcript file do not exist', logf)
      } else {
        stop('Reference transcript file do not exist')
      }
    }
    
    # use file extension format if provided by user, else try to extract file format from filename
    if (!is.null(query_format)) {
      fileformat = query_format
    } else {
      fileformat = tail(unlist(strsplit(ref, '\\.')), "1")
      if (!fileformat%in%c('gff3','gff','gtf','bed')) {
        if (makefile == FALSE) {
          stopLog('Incorrect reference file extension format', logf)
        } else {
          stop('Incorrect reference file extension format')
        }
      } 
    }
    basicGRanges = rtracklayer::import(ref, format = fileformat)
  }
  
  
  # testing and matching gene_ids
  if (makefile == FALSE) {
    infoLog('Checking and matching gene_ids...', logf, quiet)
  } else {
    message('Checking and matching gene_ids...')
  }
  
  # preserve gene_ids from input GTF
  elementMetadata(inputGRanges)$old_gene_id <- elementMetadata(inputGRanges)$gene_id
  
  # get a dataframe of all gene_ids
  #  depending on the input arguments, this dataframe is constructed differently
  #  dataframe will contain headers gene_id, new_id, match_level and might contain secondary
  #     gene_id is a reference to the orignal gene_id
  #     new_id is the newly changed id
  #     match_level describes the level at which the id has been changed
  #     0 -> ids are found in ref; 
  #     1 -> ids are matched by secondary_gene_id
  #     2 -> ids are matched by appending ENS... suffix
  #     3 -> ids are matched by secondary_gene_id followed by appending ENS... suffix
  #     4 -> ids are matched by matching overlapping coordinates
  #     5 -> id could not be matched and will be skipped from analysis
  if (!is.null(primary_gene_id) & !is.null(secondary_gene_id)) {
    unique_ids = elementMetadata(inputGRanges) %>% as.data.frame() %>%
      dplyr::select(gene_id = primary_gene_id, secondary = secondary_gene_id)
  } else if (!is.null(primary_gene_id)) {
    unique_ids = elementMetadata(inputGRanges) %>% as.data.frame() %>%
      dplyr::select(gene_id = primary_gene_id)
  } else {
    unique_ids = elementMetadata(inputGRanges) %>% as.data.frame() %>%
      dplyr::select(gene_id = gene_id)
  }
  unique_ids = unique_ids %>% 
    dplyr::mutate(new_id = gene_id, match_level = 0)
  
  # count number of non standard ID before correction
  nonstand_id_1 = nrow(elementMetadata(inputGRanges) %>% as.data.frame() %>%
                         dplyr::select(gene_id = gene_id) %>%
                         dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
                         dplyr::distinct())
  
  # There are 3 correction steps below.
  # correction 1 and 2 will initiate only if arguments are provided
  # correction 3 will always be activated, but may give false positive output
  # note that in each step, ONLY the non-standard gene_ids will be corrected
  
  # correction 1: replace primary_gene_id with secondary_gene_id, IF both args are provided
  if (!is.null(primary_gene_id) & !is.null(secondary_gene_id)) {
    if (makefile == FALSE) {
      infoLog(sprintf('-> Attempting to replace %s with %s...', primary_gene_id, secondary_gene_id),
              logf, quiet)
    } else {
      message(sprintf('-> Attempting to replace %s with %s...', primary_gene_id, secondary_gene_id))
    }
    
    # count number of non-standard ids before matching
    countsbefore = nrow(unique_ids %>%
                          dplyr::select(gene_id = new_id) %>%
                          dplyr::distinct() %>%
                          dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))

    # replace nonstandard primary_gene_ids from inputGRanges with its secondary_gene_ids
    unique_ids = unique_ids %>% 
      dplyr::mutate_at(vars(secondary), funs(ifelse(new_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE, ., NA))) %>%
      dplyr::mutate_at(vars(match_level), funs(ifelse(!is.na(secondary), 1, .))) %>%
      dplyr::mutate_at(vars(new_id), funs(ifelse(!is.na(secondary), as.character(secondary), .))) %>%
      dplyr::select(gene_id, new_id, match_level)
    
    # count number of non-standard ids after matching
    countsafter = nrow(unique_ids %>%
                         dplyr::select(gene_id = new_id) %>%
                         dplyr::distinct() %>%
                         dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))
    
    # report number of IDs corrected
    if (makefile == FALSE) {
      infoLog(sprintf('-> %s gene IDs matched', (countsbefore - countsafter)), logf, quiet)
    } else {
      message(sprintf('-> %s gene IDs matched', (countsbefore - countsafter)))
    }
  }
  
  # correction 2: replace primary_gene_id with basic gene ID IF:
  # at least primary_gene_id is provided and if it starts with 'ENS'
  if (!is.null(primary_gene_id)) {
    if (makefile == FALSE) {
      infoLog('--> Attempting to match ensembl gene_ids...', logf, quiet)
    } else {
      message('--> Attempting to match ensembl gene_ids...')
    }
    
    # count number of non-standard ids after matching
    countsbefore = nrow(unique_ids %>%
                          dplyr::select(gene_id = new_id) %>%
                          dplyr::distinct() %>%
                          dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))
                             
                              
    # make dataframe of non_standard primary_gene_ids starting with "ENS"
    gene_id_df = unique_ids %>%
      dplyr::select(gene_id = new_id) %>%
      dplyr::distinct() %>%
      dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
      dplyr::filter(startsWith(gene_id, "ENS") == TRUE)
    
    # proceed with correction if there is an "ENS" gene_id annotation
    if (nrow(gene_id_df) > 0) {
      # make dataframe of annotated gene_id starting with "ENS"
      basic_gene_id_df = elementMetadata(basicGRanges) %>% as.data.frame() %>%
        dplyr::select(ens_id = gene_id) %>%
        dplyr::distinct() %>%
        dplyr::filter(startsWith(ens_id, "ENS") == TRUE)
      
      # append suffix of ENS gene_ids on both dataframes 
      # eg. ENSMUST00000029780.11 to ENSMUST00000029780
      gene_id_df$appended_id <-  sapply(strsplit(gene_id_df$gene_id, "\\."), `[`, 1)
      basic_gene_id_df$appended_id <-  sapply(strsplit(basic_gene_id_df$ens_id, "\\."), `[`, 1)
      
      # combine both data frames based on common 'appended_id' and filter lines with empty 'ens_id'
      gene_id_df = dplyr::left_join(gene_id_df, basic_gene_id_df, by = 'appended_id') %>%
        dplyr::filter(!is.na(ens_id))
      
      # replace new_id on unique_ids dataframe with ens_id, if match is found
      unique_ids = unique_ids %>% left_join(gene_id_df, by=c('new_id'='gene_id')) %>%
        dplyr::mutate_at(vars(match_level), funs(ifelse(!is.na(ens_id), as.integer(.)+1, .))) %>%
        dplyr::mutate_at(vars(new_id), funs(ifelse(!is.na(ens_id), as.character(ens_id), .))) %>%
        dplyr::select(gene_id, new_id, match_level)
      
      # count number of non-standard ids after matching
      countsafter = nrow(unique_ids %>%
                           dplyr::select(gene_id = new_id) %>%
                           dplyr::distinct() %>%
                           dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))
      
      # report number of IDs corrected
      if (makefile == FALSE) {
        infoLog(sprintf('--> %s gene IDs matched', (countsbefore - countsafter)), logf, quiet)
      } else {
        message(sprintf('--> %s gene IDs matched', (countsbefore - countsafter)))
      }
      
    } else {
      # check if there is any ENS ids to begin with
      ENSids = unique_ids %>%
        dplyr::distinct() %>%
        dplyr::filter(startsWith(gene_id, "ENS") == TRUE)
      
      if (nrow(ENSids) > 0) {
        if (makefile == FALSE) {
          infoLog('--> All ensembl gene ids have been matched', logf, quiet)
        } else {
          message('--> All ensembl gene ids have been matched')
        }
      } else {
        if (makefile == FALSE) {
          warnLog('--> No ensembl gene ids found in query', logf, quiet)
        } else {
          message('--> No ensembl gene ids found in query')
        }
      }
    }
  }
  
  # correction 3: correct gene_ids by finding overlapping regions.
  if (makefile == FALSE) {
    infoLog('---> Attempting to correct gene_ids by finding overlapping coordinates...', logf, quiet)
  } else {
    message('---> Attempting to correct gene_ids by finding overlapping coordinates...')
  }
  
  # count number of unmatched ids before matching
  countsbefore = nrow(unique_ids %>%
                        dplyr::select(gene_id = new_id) %>%
                        dplyr::distinct() %>%
                        dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))
  
  # make dataframe of non_standard gene_ids
  gene_id_df = unique_ids %>%
    dplyr::select(gene_id = new_id, test_id = gene_id) %>%
    dplyr::distinct() %>%
    dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
    as.data.frame()
  
  # this is the main worker of this script
  #   function will retrieve first gene_id from gencode that overlaps with first exon of query
  #   the below parameters give a true positive rate of 81% 
  if (nrow(gene_id_df) > 0) {
    
    group <- rep(1:clusters, length.out = nrow(gene_id_df))
    gene_id_df <- bind_cols(tibble(group), gene_id_df)
    cluster <- create_cluster(cores = clusters, quiet = TRUE)
    
    parallel_df = gene_id_df %>% partition(group, cluster = cluster)
    parallel_df %>%
      # Assign libraries
      cluster_library("GenomicRanges") %>%
      cluster_library("dplyr") %>%
      # Assign values (use this to load functions or data to each core)
      cluster_assign_value("inputGRanges", inputGRanges) %>%
      cluster_assign_value("basicGRanges", basicGRanges)
    
    gene_id_df <- parallel_df %>% # Use by_group party_df
      do(dplyr::mutate(., ref_gene_id = sapply(.$test_id, function(x) {
                      
                      getfirstOverlap = findOverlaps(
                        inputGRanges[mcols(inputGRanges)$gene_id == x &
                                       mcols(inputGRanges)$type == "exon"][1],
                        basicGRanges, 
                        select = "first")
                      
                      return(ifelse(!is.na(getfirstOverlap), 
                                    mcols(basicGRanges[getfirstOverlap])$gene_id,
                                    NA))
                    }))) %>%
      collect() %>% # Special collect() function to recombine partitions
      as.data.frame() %>%
      dplyr::filter(!is.na(ref_gene_id))
    

    
    # replace gene_id on assembled transcript with predicted gene_id from overlap
    unique_ids = unique_ids %>% left_join(gene_id_df, by=c('new_id'='gene_id')) %>%
      dplyr::mutate_at(vars(match_level), funs(ifelse(!is.na(ref_gene_id), 4, .))) %>%
      dplyr::mutate_at(vars(new_id), funs(ifelse(!is.na(ref_gene_id), as.character(ref_gene_id), .))) %>%
      dplyr::select(gene_id, new_id, match_level)
    
    # count number of unmatched ids after matching
    countsafter = nrow(unique_ids %>%
                         dplyr::select(gene_id = new_id) %>%
                         dplyr::distinct() %>%
                         dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))
    
    # report number of IDs corrected
    if (makefile == FALSE) {
      infoLog(sprintf('---> %s gene IDs matched', (countsbefore - countsafter)), logf, quiet)
    } else {
      message(sprintf('---> %s gene IDs matched', (countsbefore - countsafter)))
    }
  } else {
    if (makefile == FALSE) {
      infoLog(sprintf('---> Skipped, all gene IDs matched', (countsbefore - countsafter)), logf, quiet)
    } else {
      message(sprintf('---> Skipped, all gene IDs matched', (countsbefore - countsafter)))
    }
  }
  

  # count remaining number of non_standard gene_ids and number of corrected ids
  nonstand_id_list = unique_ids %>%
    dplyr::select(gene_id = new_id) %>%
    dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
    dplyr::distinct() %>%
    dplyr::mutate(status = TRUE)
  nonstand_id_2 = nrow(nonstand_id_list)
  corrected_ids = nonstand_id_1 - nonstand_id_2
  
  if (nonstand_id_2 > 0) {
    
    # if there are remaining unmatched gene_ids, update match_level to 5.
    #   these transcripts will be skipped from further analysis
    unique_ids = unique_ids %>% left_join(nonstand_id_list, by=c('new_id'='gene_id')) %>%
      dplyr::mutate_at(vars(match_level), funs(ifelse(!is.na(status), 5, .))) %>%
      dplyr::select(gene_id, new_id, match_level)
  }

  # update gene_id and add match_level to inputGRanges
  mcols(inputGRanges)$gene_id = unique_ids$new_id
  mcols(inputGRanges)$match_level = unique_ids$match_level

  # report pre-testing analysis and return inputGRanges
  if (makefile == FALSE) {
    infoLog(sprintf('Total gene_ids corrected: %s', corrected_ids), logf, quiet)
    infoLog(sprintf('Remaining number of non-standard gene_ids: %s', nonstand_id_2), logf, quiet)
    if (nonstand_id_2 > 0) {
      warnLog('Transcripts with non-standard gene_ids will be skipped', logf, quiet)
    }
  } else {
    message(sprintf('Total gene_ids corrected: %s', corrected_ids))
    message(sprintf('Remaining number of non-standard gene_ids: %s', nonstand_id_2))
    if (nonstand_id_2 > 0) {
      message('Transcripts with non-standard gene_ids will be skipped')
    }
  }
  if (makefile == TRUE) {
    rtracklayer::export(inputGRanges, outputfile)
    message(sprintf('Done. GTF saved as %s in current working directory', outputfile))
  } else {
    return(inputGRanges)
  }
}

