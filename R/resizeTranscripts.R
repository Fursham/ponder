resizeTranscripts <- function(x, headlength = 0, taillength = 0) {
  
  if (sum(width(x)) < (headlength + taillength)) {
    stop('Appending length is larger than size of transcript')
  }
  
  # retrieve important information
  strand = as.character(strand(x))[1]
  x = x %>% as.data.frame() %>%
    dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start)) %>%
    dplyr::mutate(fwdcumsum = as.integer(cumsum(width) - headlength), 
                  revcumsum = rev(cumsum(rev(width))) - taillength) %>%
    dplyr::filter(fwdcumsum > 0 & revcumsum > 0)
  
  if(strand == '-'){
    x[1,]$end = x[1,]$start + x[1,]$fwdcumsum -1
    x[nrow(x),]$start = x[nrow(x),]$end - x[nrow(x),]$revcumsum + 1
  } else{
    x[1,]$start = x[1,]$end - x[1,]$fwdcumsum + 1
    x[nrow(x),]$end = x[nrow(x),]$start + x[nrow(x),]$revcumsum - 1
  }
  x = x %>% dplyr::select(-fwdcumsum,-revcumsum) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  return(x)
}