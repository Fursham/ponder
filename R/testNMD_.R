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
#' 
#' 
testNMD_ <- function(queryTranscript, queryCDS, distance_stop_EJ = 50, other_features = FALSE, fasta){
  
  # prepare output list
  output = list(is_NMD = as.logical(FALSE), 
                dist_to_lastEJ = as.numeric(0),
                num_of_down_EJs = as.numeric(0),
                dist_to_downEJs = as.numeric(0))
  if (other_features == TRUE) {
    output = modifyList(output, list(threeUTRlength = as.numeric(NA),
                                     uORF = as.character(NA), 
                                     uATG = as.character(NA), 
                                     uATG_frame = as.character(NA)))
  }
  

  # sort queryCDS, by exon order (just in case)
  strand = as.character(BiocGenerics::strand(queryCDS))[1]
  queryCDS = BiocGenerics::sort(queryCDS, decreasing = strand == '-')
  queryTranscript = BiocGenerics::sort(queryTranscript, decreasing = strand == '-')
  
  # test if query is NMD sensitive
  #   get distance of last EJC from start of transcript
  #   disjoin will create a new GRanges that will separate the queryTranscript
  #   into discernable 5UTR, ORF and 3UTR
  #   we can then try to use the new GRanges to infer NMD susceptibility
  lengthtolastEJ = sum(head(BiocGenerics::width(queryTranscript),-1)) 
  disjoint = BiocGenerics::append(queryCDS,queryTranscript) %>%
    GenomicRanges::disjoin(with.revmap = T) %>%
    BiocGenerics::sort(decreasing = strand == '-') %>%
    as.data.frame() %>%
    dplyr::mutate(disttolastEJ = lengthtolastEJ - cumsum(width))
  
  # retrieve index of last ORF segment and determine if there are 
  # more than 1 exons after stop codon
  stopcodonindex = max(which(lengths(disjoint$revmap) == 2))
  output$dist_to_lastEJ = disjoint[stopcodonindex,]$disttolastEJ
 
  # report number of downstream EJ and its distance to PTC
  output$num_of_down_EJs = nrow(disjoint) - stopcodonindex -1
  downEJCdf = disjoint %>% dplyr::filter(dplyr::row_number() >= stopcodonindex+1) %>%
    dplyr::filter(disttolastEJ >= 0) %>%
    dplyr::mutate(disttoPTC = cumsum(width))
  output$dist_to_downEJs = paste(downEJCdf$disttoPTC, collapse = ',')
  
  if(output$dist_to_lastEJ > distance_stop_EJ){
    #annotated transcript as NMD if dist_to_lastEJ is NMD triggering
    output$is_NMD = TRUE  
  }
  

  # test for other features if activated
  if (other_features == TRUE) {
    if (missing(fasta)) {
      stop('Please provide fasta sequence')
    }
    
    # report number of downstream EJ and its distance to PTC
    output$num_of_down_EJs = nrow(disjoint) - stopcodonindex -1
    downEJCdf = disjoint %>% dplyr::filter(dplyr::row_number() >= stopcodonindex+1) %>%
      dplyr::filter(disttolastEJ >= 0) %>%
      dplyr::mutate(disttoPTC = cumsum(width))
    output$dist_to_downEJs = paste(downEJCdf$disttoPTC, collapse = ',')
    
    # obtain position of start/stop codons in disjointed GRanges
    startcodonindex = min(which(lengths(disjoint$revmap) == 2))
    stopcodonindex = max(which(lengths(disjoint$revmap) == 2))
    
    if(startcodonindex >1){
      fiveUTRindex = startcodonindex -1
    } else{
      fiveUTRindex = 1
    }
    

    if(stopcodonindex < nrow(disjoint)){
      threeUTRindex = stopcodonindex + 1
    } else {
      threeUTRindex = stopcodonindex 
    }
    
    # obtain cumsumlength of disjointed GRanges and obtain UTR lengths
    disjoint = disjoint %>% as.data.frame() %>%
      dplyr::mutate(tmp.fwdcumsum = cumsum(width),
                    tmp.revcumsum = rev(cumsum(rev(width))))
    fiveUTRlength = disjoint[fiveUTRindex,]$tmp.fwdcumsum
    output$threeUTRlength = disjoint[threeUTRindex,]$tmp.revcumsum

    # test for uORF on 5'UTR
    # get sequence of 5'UTR
    fiveUTRGRanges = resizeGRangesTranscripts(queryTranscript, end = sum(BiocGenerics::width(queryTranscript)) - fiveUTRlength)
    list_startstopcodons = Biostrings::DNAStringSet(c("ATG","TAA", "TAG", "TGA"))
    pdict_startstopcodons = Biostrings::PDict(list_startstopcodons)
    
    # this part will test the presence of uORFs and uATGs in the 5UTR
    fiveUTRseq = unlist(BSgenome::getSeq(fasta, fiveUTRGRanges))
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
    gr = uORFuATG %>% dplyr::mutate(seqnames = 1) %>% GenomicRanges::makeGRangesFromDataFrame()
    nonOverlapsuORFuATG = uORFuATG[BiocGenerics::unique(GenomicRanges::findOverlaps(gr, type = "any", select = "first")),] %>%
      dplyr::slice(1:ifelse('uATG'%in%group_name, min(which('uATG' == group_name)), dplyr::n()))

    # retrieve GRanges for the uORFs and uATGs above
    uORFGranges = do.call('c', base::mapply(function(x,y,z,a){
      start = x - 1
      end = length(fiveUTRseq) - y
      startcodoninGRanges = range(resizeGRangesTranscripts(fiveUTRGRanges, start, end))
      S4Vectors::mcols(startcodoninGRanges)$group_name = z
      S4Vectors::mcols(startcodoninGRanges)$frame = a
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
