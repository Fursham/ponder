#' GFF3 to GTF converter
#'
#' @param gffGRanges Input GFF3 GRanges object
#'
#' @return Outputs GTF GRanges object with mandatory gene_id and transcript_id headers
#' @export
#'
gff3togtfconvert <- function(gffGRanges) {
  
  if(!'Parent' %in% names(mcols(gffGRanges))){
    stopLog('Please check that GFF3 file contain Parent information')
  }
  
  gffGRanges = gffGRanges %>% as.data.frame() %>%
    dplyr::mutate(gene_id = ifelse(type == 'gene', ID, NA)) %>% 
    dplyr::mutate(transcript_id = ifelse(type == 'transcript', ID, NA)) %>%
    dplyr::mutate(Parent = as.character(Parent))
    
    
  gtfGRanges = gffGRanges %>% as.data.frame() %>%
    dplyr::mutate(gene_id = ifelse(Parent %in% gffGRanges$gene_id, Parent, gene_id)) %>% 
    dplyr::mutate(transcript_id = ifelse(Parent %in% gffGRanges$transcript_id, Parent, transcript_id)) %>% 
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(desc(width)) %>%
    tidyr::fill() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(start) %>%
    dplyr::select(-Parent) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  return(gtfGRanges)
}