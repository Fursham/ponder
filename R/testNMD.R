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
  output = list(is_NMD = as.logical(FALSE), dist_to_lastEJ = as.numeric(0))
  if (other_features == TRUE) {
    output = modifyList(output, list(uORF = as.character(NA), 
                                     threeUTR = as.numeric(NA), 
                                     uATG = as.character(NA), 
                                     uATG_frame = as.character(NA)))
  }
  

  # sort queryCDS, by exon order (just in case)
  strand = as.character(strand(queryCDS))[1]
  queryCDS = sort(queryCDS, decreasing = strand == '-')
  queryTranscript = sort(queryTranscript, decreasing = strand == '-')
  
  disjoint = BiocGenerics::append(queryCDS,queryTranscript) %>%
    GenomicRanges::disjoin(with.revmap = T) %>%
    sort(decreasing = strand == '-')
  stopcodonindex = max(which(lengths(mcols(disjoint)$revmap) == 2))
  num_of_exons_aft_stop = length(disjoint) - stopcodonindex
  if(num_of_exons_aft_stop > 1){
    output$dist_to_lastEJ = width(disjoint[stopcodonindex+1])
    if(output$dist_to_lastEJ > distance_stop_EJ){
      output$is_NMD = TRUE
    }
  }



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
