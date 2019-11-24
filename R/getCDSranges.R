getCDSranges <- function(query, fiveUTRlength, threeUTRlength, gene_id, transcript_id){
  
  CDSranges = resizeTranscripts(query, fiveUTRlength, threeUTRlength)
  CDSranges = CDSranges %>% as.data.frame() %>%
    dplyr::mutate(type = 'CDS', gene_id = gene_id, transcript_id = transcript_id) %>%
    dplyr::mutate(phase = cumsum(width%%3)%%3)
  CDSranges$phase = c(0, head(augmentedCDS$phase, - 1))
  CDSranges = makeGRangesFromDataFrame(augmentedCDS, keep.extra.columns = TRUE)
  
  return(list(ORF_considered = CDSranges))
}