testAS <- function(tx1, tx2){
  
  #prepare output
  output = as.character(NA)
  
  if(!any(tx1 %in% tx2)){
    return()
  }
  
  # get information on transcripts
  tx1index = c(1:length(tx1))
  tx2index = c(length(tx1)+length(tx1)+length(tx2))
  tx1ends = list('start' = 1, 
                 'end' = length(tx1))
  tx2ends = list('start' = length(tx1)+1, 
                 'end' = length(tx1)+length(tx2))
  txcomb = c(tx1ends, tx2ends)
  strand = strand(tx1)[1] %>% as.character()

  # combine transcripts and disjoin
  disjoint = BiocGenerics::append(tx1,tx2) %>%
    GenomicRanges::disjoin(with.revmap = T) %>%
    as.data.frame() %>%
    dplyr::mutate(type = ifelse(lengths(revmap)==2, 'cons','alt')) 
  disjoint[1:min(which(lengths(disjoint$revmap)==2))-1,]$type = 'up'
  disjoint[(max(which(lengths(disjoint$revmap)==2))+1):nrow(disjoint),]$type = 'down'
  
  disjoint = disjoint %>%
    dplyr::mutate(source = ifelse(type != 'cons' & revmap%in%tx1index, 1, 0)) %>%
    dplyr::mutate(source = ifelse(type != 'cons' & revmap%in%tx2index, 2, source)) %>%
    dplyr::mutate(typesource = paste0(type,'_',source)) %>%
    dplyr::mutate(typelag = dplyr::lag(type, default = 'start'), typelead = dplyr::lead(type, default = 'end')) %>%
    dplyr::mutate(upcoord = dplyr::lag(end, default = .$start[1]), downcoord = dplyr::lead(start, default = .$end[dplyr::n()])) %>%
    dplyr::mutate(updiff = start - upcoord, downdiff = downcoord-end)
    
    
    dplyr::mutate(type = ifelse(type == 'alt' & revmap%in%txcomb, 'ends', type)) %>%
    dplyr::mutate(type = ifelse(type == 'alt', 'internal', type)) %>%
    dplyr::mutate(typesource = paste0(type,'_',source)) %>%
    dplyr::mutate(sourcelag = dplyr::lag(typesource, default = 'start'), sourcelead = dplyr::lead(typesource, default = 'end')) %>%
    dplyr::mutate(upcoord = dplyr::lag(end, default = .$start[1]), downcoord = dplyr::lead(start, default = .$end[dplyr::n()])) %>%
    dplyr::mutate(updiff = start - upcoord, downdiff = downcoord-end)
  
  disjoint = disjoint %>%
    dplyr::filter(type != 'cons') %>%
    dplyr::mutate(AS = NA) %>%
    dplyr::mutate(AS = ifelse(sourcelag == 'start' & sourcelead == 'cons_0', 'ATS',AS)) %>%
    dplyr::mutate(AS = ifelse(sourcelag == 'start' & sourcelead == paste0('alt_', ifelse(source == 1,2,1)), 'AF', AS)) %>%
    dplyr::mutate(AS = ifelse(sourcelag != 'cons_0' & sourcelead == 'cons_0', 'AF',AS)) %>%
    
    dplyr::mutate(AS = ifelse(type == 'ends' & sourcelag == 'cons_0' & sourcelead == 'cons_0' & downdiff == 0, 'APA',AS)) %>%
    dplyr::mutate(AS = ifelse(type == 'ends' & (sourcelag != 'cons_0' | sourcelead != 'cons_0') & updiff == 0, 'ATS',AS)) %>%
    dplyr::mutate(AS = ifelse(type == 'ends' & (sourcelag != 'cons_0' | sourcelead != 'cons_0') & downdiff == 0, 'APA',AS)) %>%

    
}

