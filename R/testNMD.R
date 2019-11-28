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
    output = modifyList(output, list(threeUTRlength = as.numeric(NA),
                                     uORF = as.character(NA), 
                                     uATG = as.character(NA), 
                                     uATG_frame = as.character(NA)))
  }
  

  # sort queryCDS, by exon order (just in case)
  strand = as.character(strand(queryCDS))[1]
  queryCDS = sort(queryCDS, decreasing = strand == '-')
  queryTranscript = sort(queryTranscript, decreasing = strand == '-')
  
  # test if query is NMD sensitive
  #   disjoin will create a new GRanges that will separate the queryTranscript
  #   into discernable 5UTR, ORF and 3UTR
  #   we can then try to use the new GRanges to infer NMD susceptibility
  disjoint = BiocGenerics::append(queryCDS,queryTranscript) %>%
    GenomicRanges::disjoin(with.revmap = T) %>%
    sort(decreasing = strand == '-')
  
  # retrieve index of last ORF segment and determine if there are 
  # more than 1 exons after stop codon
  stopcodonindex = max(which(lengths(mcols(disjoint)$revmap) == 2))
  num_of_exons_aft_stop = length(disjoint) - stopcodonindex
  
  # if more than 1 segment is found, test for distance to last EJ
  if(num_of_exons_aft_stop > 1){
    output$dist_to_lastEJ = width(disjoint[stopcodonindex+1])
    if(output$dist_to_lastEJ > distance_stop_EJ){
      #annotated transcript as NMD if dist_to_lastEJ is NMD triggering
      output$is_NMD = TRUE  
    }
  }



  # test for other features if activated
  if (other_features == TRUE) {
    if (missing(fasta)) {
      stop('Please provide fasta sequence')
    }
    
    # obtain position of start/stop codons in disjointed GRanges
    startcodonindex = min(which(lengths(mcols(disjoint)$revmap) == 2))
    stopcodonindex = max(which(lengths(mcols(disjoint)$revmap) == 2))
    
    # obtain cumsumlength of disjointed GRanges and obtain UTR lengths
    disjoint = disjoint %>% as.data.frame() %>%
      dplyr::mutate(tmp.fwdcumsum = cumsum(width),
                    tmp.revcumsum = rev(cumsum(rev(width))))
    fiveUTRlength = disjoint[startcodonindex-1,]$tmp.fwdcumsum
    output$threeUTRlength = disjoint[stopcodonindex+1,]$tmp.revcumsum

    # test for uORF on 5'UTR
    # get sequence of 5'UTR
    fiveUTRGRanges = resizeGRangesTranscripts(queryTranscript, end = sum(width(queryTranscript)) - fiveUTRlength)
    list_startstopcodons = Biostrings::DNAStringSet(c("ATG","TAA", "TAG", "TGA"))
    pdict_startstopcodons = Biostrings::PDict(list_startstopcodons)
    
    # this part will test the presence of uORFs and uATGs in the 5UTR
    fiveUTRseq = unlist(Biostrings::getSeq(fasta, fiveUTRGRanges))
    allmatches = Biostrings::matchPDict(pdict_startstopcodons, fiveUTRseq) %>%
      as.data.frame()
    
    # return if 5UTR contain no start/stop codons
    if(nrow(allmatches) == 0){
      return(output)
    }
    
    # this code will generate a dataframe of start/stop coordinates of
    # uORFs and uATGs
    uORFuATG = allmatches %>%
      dplyr::mutate(group_name = ifelse(group == 1, 'start', 'stop'),
                    frame = (length(fiveUTRseq) - end)%%3) %>%
      dplyr::arrange(start) %>%
      dplyr::group_by(frame) %>%
      dplyr::mutate(shiftype = dplyr::lag(group_name, default = 'stop')) %>%
      dplyr::filter(group_name != shiftype) %>%
      dplyr::mutate(group_name = ifelse(dplyr::n()%%2 != 0 & dplyr::row_number()==dplyr::n(),
                           'uATG', group_name)) %>%
      dplyr::mutate(end = ifelse(group_name == 'start', NA, end)) %>%
      tidyr::fill(end, .direction = 'up') %>%
      dplyr::filter(group_name != 'stop') %>%
      dplyr::mutate(group_name = ifelse(group_name == 'start', 'uORF', group_name),
                    width = end - start + 1) %>% 
      dplyr::ungroup() %>%
      dplyr::arrange(start) %>%
      dplyr::select(-shiftype)
    
    # return if no uORFs or uATGs are found
    if(nrow(uORFuATG) == 0){
      return(output)
    }
    
    # this code will return non-overlapping uORFs and stops at the first uATG
    gr = uORFuATG %>% dplyr::mutate(seqnames = 1) %>% makeGRangesFromDataFrame()
    nonOverlapsuORFuATG = uORFuATG[unique(findOverlaps(gr, type = "any", select = "first")),] %>%
      dplyr::slice(1:ifelse('uATG'%in%group_name, min(which('uATG' == group_name)), dplyr::n()))

    # retrieve GRanges for the uORFs and uATGs above
    uORFGranges = do.call('c', base::mapply(function(x,y,z,a){
      start = x - 1
      end = length(fiveUTRseq) - y
      startcodoninGRanges = range(resizeGRangesTranscripts(fiveUTRGRanges, start, end))
      mcols(startcodoninGRanges)$group_name = z
      mcols(startcodoninGRanges)$frame = a
      return(startcodoninGRanges)
    }, nonOverlapsuORFuATG$start, nonOverlapsuORFuATG$end, 
    nonOverlapsuORFuATG$group_name, nonOverlapsuORFuATG$frame)) %>% 
      as.data.frame()
    
    if('uORF'%in%uORFGranges$group_name){
      uORF = uORFGranges %>% 
        dplyr::filter(group_name == 'uORF') %>%
        dplyr::mutate(coord = paste0(start, '-', end)) %>%
        dplyr::select(coord)
      output$uORF = paste(uORF$coord, collapse = '|')
    }
    
    if('uATG'%in%uORFGranges$group_name){
      uATG = uORFGranges %>% 
        dplyr::filter(group_name == 'uATG') %>%
        dplyr::mutate(coord = paste0(start, '-', end)) %>%
        dplyr::select(coord, frame)
      output$uATG = paste(uATG$coord, collapse = '|')
      output$uATG_frame = paste(uATG$frame, collapse = '|')
    }
  }
  # return output
  return(output)
}
