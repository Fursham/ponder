searchuORFs <- function(exons, cds, fasta,
                           ORFlength = 21, which = NULL){
  
  # catch missing args
  mandargs <- c('exons', 'cds','fasta')
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    stop(paste("missing values for", 
               paste(setdiff(mandargs, passed), collapse=", ")))
  }
  
  # check if exons and cds are GRlist
  if(!is(exons,'GRangesList') | !is(cds,'GRangesList')){
    txtype = is(exons)[1]
    cdstype = is(cds)[1]
    
    error = c(!is(exons,'GRangesList') , !is(cds,'GRangesList'))
    obj = paste(c('exons','cds')[error], collapse = ',')
    types = paste(c(is(exons)[1], is(cds)[1])[error], collapse = ',')

    stop(sprintf('Incompatile input object. %s is type %s respectively',
                 obj, types))
  }
  
  # catch unmatched seqlevels
  if(GenomeInfoDb::seqlevelsStyle(exons) != GenomeInfoDb::seqlevelsStyle(cds)){
    stop('exons and cds has unmatched seqlevel styles. try matching using matchSeqLevels function')
  }
  
  # test for uORF and uATG
  totest = names(exons)
  # check for missing cds and return warnings/errors
  totest = totest[totest %in% names(cds)]
  if(length(totest)==0){
    stop('all tx have missing cds info. please ensure tx and cds names match')
  } 
  if(length(totest) < length(exons)){
    skiptest = length(exons) - length(totest)
    rlang::warn(sprintf('%s tx(s) have missing cds info and have been skipped',
                        skiptest))
  }
  out = BiocParallel::bplapply(totest,  function(x){
    report = getuORFuATG_(exons[x], cds[x], fasta, ORFlength)
    return(report %>% as.data.frame())
  }, BPPARAM = BiocParallel::MulticoreParam()) %>%
    dplyr::bind_rows()
  
  newtx <- out %>% dplyr::select(dplyr::starts_with('tx')) %>%
    dplyr::distinct() %>%
    dplyr::rename_all(function(x){return(substr(x, start = 4, stop = nchar(x)))}) %>%
    dplyr::select(-group) %>%
    GenomicRanges::makeGRangesListFromDataFrame(split.field = 'group_name', 
                                                keep.extra.columns = T)
  
  newcds <- out %>% dplyr::select(dplyr::starts_with('cds')) %>%
    dplyr::distinct() %>%
    dplyr::rename_all(function(x){return(substr(x, start = 5, stop = nchar(x)))}) %>%
    dplyr::select(-group) %>%
    GenomicRanges::makeGRangesListFromDataFrame(split.field = 'group_name', 
                                                keep.extra.columns = T)
  
  return(list('exons' = newtx, 'cds' = newcds))
}



getuORFuATG_ <- function(txlist, cdslist, fasta, size){
  
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
    return() #return if start tx and cds starts with start codon
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
    return()
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
    dplyr::filter(width >= size) %>%
    dplyr::mutate(group_name = ifelse(group_name == 'uORF',
                                      paste0(group_name, '_',dplyr::row_number()),
                                      group_name)) %>% 
    dplyr::select(-shiftype)
  
  # return if no uORFs or uATGs are found
  if(nrow(uORFuATG) == 0){
    return()
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
    S4Vectors::mcols(thisGR)$type = 'CDS'
    S4Vectors::mcols(thisGR)$frame = a
    S4Vectors::mcols(thisGR)$phase <- data.table::shift(cumsum(width(thisGR)%%3)%%3, fill = 0)
    S4Vectors::mcols(thisGR)$newname= paste0(z, '_', names(txlist))
    return(thisGR)
  }, nonOverlapsuORFuATG$start, nonOverlapsuORFuATG$end, 
  nonOverlapsuORFuATG$group_name, nonOverlapsuORFuATG$frame)) %>% 
    as.data.frame()
  uORFgr = uORFs %>%
    makeGRangesListFromDataFrame(split.field = 'newname', 
                                 keep.extra.columns = T)
  
  # made new tx GRangesList
  txlistnew = rep(txlist,length(uORFgr))
  names(txlistnew) = names(uORFgr)
  
  # get CDS with newly found uATG
  if(any(startsWith(names(uORFgr), 'uATG'))){
    uATGindex = which(startsWith(names(uORFgr),'uATG'))
    uATGgr = uORFgr[uATGindex]
    uATGtx = txlistnew[names(uATGgr)]
    uATGCDS = getCDS_(uATGtx, uATGgr, fasta)
    if(!is.null(uATGCDS)){
      uATGCDS = uATGCDS %>%
        dplyr::select(-transcript_id) %>%
        makeGRangesFromDataFrame(keep.extra.columns = T)
      S4Vectors::mcols(uATGCDS)$frame= mcols(uATGgr[[1]])$frame
      uORFgr[[uATGindex]] = uATGCDS
    }
  }
  return(list('tx' = txlistnew, 'cds' = uORFgr))
}
