uORFuATG_ <- function(txlist, cdslist, fasta){
  
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
  
  #extract GR from GRL
  tx = txlist[[1]]
  cds = cdslist[[1]]
  
  # obtain strand and sort, just in case
  strand = as.character(BiocGenerics::strand(cds))[1]
  txwidth = sum(BiocGenerics::width(tx))
  cds = BiocGenerics::sort(cds, decreasing = strand == '-')
  tx = BiocGenerics::sort(tx, decreasing = strand == '-')
  
  # test if query is NMD sensitive
  #   get distance of last EJC from start of transcript
  #   disjoin will create a new GRanges that will separate the tx
  #   into discernable 5UTR, ORF and 3UTR
  #   we can then try to use the new GRanges to infer NMD susceptibility
  disjoint = BiocGenerics::append(cds,tx) %>%
    GenomicRanges::disjoin(with.revmap = T) %>%
    BiocGenerics::sort(decreasing = strand == '-') %>%
    as.data.frame() %>%
    dplyr::mutate(tmp.fwdcumsum = cumsum(width),
                  tmp.revcumsum = rev(cumsum(rev(width))))
  

  # obtain position of start codons in disjointed GRanges and fiveUTR
  startcodonindex = min(which(lengths(disjoint$revmap) == 2))

  if(startcodonindex >1){
    fiveUTRindex = startcodonindex -1
  } else{
    fiveUTRindex = 1
  }
  fiveUTRlength = disjoint[fiveUTRindex,]$tmp.fwdcumsum


  
  # test for uORF on 5'UTR
  # get sequence of 5'UTR
  fiveUTRGRanges = resizeGRangesTranscripts(tx, end = txwidth - fiveUTRlength)
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
    dplyr::mutate(group_name = ifelse(group_name == 'uORF',
                                      paste0(dplyr::row_number(),group_name),
                                      group_name)) %>% 
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
  uORFs = do.call('c', base::mapply(function(x,y,z,a){
    start = x - 1
    end = length(fiveUTRseq) - y
    thisGR = resizeGRangesTranscripts(fiveUTRGRanges, start, end) 
    S4Vectors::mcols(thisGR)$group_name = z
    S4Vectors::mcols(thisGR)$frame = a
    S4Vectors::mcols(thisGR)$newname= paste0(z, '_', names(txlist))
    return(thisGR)
  }, nonOverlapsuORFuATG$start, nonOverlapsuORFuATG$end, 
  nonOverlapsuORFuATG$group_name, nonOverlapsuORFuATG$frame)) %>% 
    as.data.frame()
  uORFgr = uORFs %>%
    makeGRangesListFromDataFrame(split.field = 'newname', 
                                 keep.extra.columns = T)
  
  
  if('uATG'%in%uORFs$group_name){
    uATGgr = uORFgr[[which('uATG' == uORFs$group_name)]]
    disjoint = BiocGenerics::append(uATGgr,tx) %>%
      GenomicRanges::disjoin(with.revmap = T) %>%
      BiocGenerics::sort(decreasing = strand == '-')
  }
  
  
  ############
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
