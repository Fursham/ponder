#' Classify alternative splicing between transcripts
#'
#' @param tx1 GRanges object of exons by transcript 1
#' @param tx2 GRanges object of exons by transcript 2
#'
#' @return GRanges object of alternative segments
#' @export
#'
#' @examples
classifyAS <- function(tx1, tx2){
  
  # input checks
  if (missing(tx1) | missing(tx2)) {
    stop('Please provide input GRanges objects')
  } else if(!any(tx1 %over% tx2)){
    stop('Transcripts have no overlap')
  }
  
  # get information on transcripts
  tx1index = c(1:length(tx1))
  tx2index = c(length(tx1)+length(tx1)+length(tx2))
  strand = strand(tx1)[1] %>% as.character()

  # combine transcripts and disjoin, annotate position of alt segments
  disjoint = BiocGenerics::append(tx1,tx2) %>%
    GenomicRanges::disjoin(with.revmap = T) %>%
    as.data.frame() %>%
    dplyr::mutate(type = ifelse(lengths(revmap)==2, 'cons','alt')) 
  disjoint[1:min(which(disjoint$type == 'cons'))-1,]$type = 'up'
  disjoint[(max(which(disjoint$type == 'cons'))+1):nrow(disjoint),]$type = 'down'
  
  # prepare dataframe for classifcation
  disjoint = disjoint %>%
    dplyr::mutate(source = ifelse(type != 'cons' & revmap%in%tx1index, 1, 0)) %>%
    dplyr::mutate(source = ifelse(type != 'cons' & revmap%in%tx2index, 2, source)) %>%
    dplyr::mutate(countersource = chartr("12", "21", source)) %>%
    dplyr::mutate(sourcedown = dplyr::lag(source, default = 0), sourceup = dplyr::lead(source, default = 0)) %>%
    dplyr::mutate(upcoord = dplyr::lag(end, default = .$start[1]), downcoord = dplyr::lead(start, default = .$end[dplyr::n()])) %>%
    dplyr::mutate(updiff = start - upcoord, downdiff = downcoord-end)
  
  # classify AS based on features
  disjoint = disjoint %>%
    dplyr::filter(type != 'cons') %>%
    dplyr::mutate(AS = NA) %>%
    dplyr::mutate(AS = ifelse(type == 'up' & downdiff == 1, 'ts',AS)) %>%
    dplyr::mutate(AS = ifelse(type == 'up' & downdiff > 1, 'FE',AS)) %>%
    dplyr::mutate(AS = ifelse(type == 'down' &updiff == 1, 'pa',AS)) %>%
    dplyr::mutate(AS = ifelse(type == 'down' & updiff > 1, 'LE',AS)) %>%
    dplyr::mutate(AS = ifelse(type == 'alt' & updiff == 1 & downdiff == 1, 'RI',AS)) %>%
    dplyr::mutate(AS = ifelse(type == 'alt' & updiff == 1 & downdiff > 1, 'SD',AS)) %>%
    dplyr::mutate(AS = ifelse(type == 'alt' & updiff > 1 & downdiff == 1, 'SA',AS)) %>%
    dplyr::mutate(AS = ifelse(type == 'alt' & updiff > 1 & downdiff > 1 & (countersource %in% c(sourcedown, sourceup)), 'ME',AS)) %>%
    dplyr::mutate(AS = ifelse(type == 'alt' & updiff > 1 & downdiff > 1 & (!countersource %in% c(sourcedown, sourceup)), 'CE',AS))
  
  # convert from alternative first exon to alternative (e.g.) if Ranges is on the neg strand
  if(strand == '-'){
    disjoint = disjoint %>%
      dplyr::mutate(AS = chartr('DAFLtspa', 'ADLFpats', AS)) %>%
      dplyr::arrange(dplyr::desc(start))
  }
  
  # convert cases depending on whether segment is found in tx1 or tx2
  #return as GRanges
  disjoint = disjoint %>% 
    dplyr::mutate(AS = ifelse(source == 2, toupper(AS), tolower(AS))) %>%
    dplyr::select(seqnames:strand, AS) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  return(disjoint)

}

