getCDSranges <- function(query, fiveUTRlength, threeUTRlength){
  
  CDSranges = resizeTranscripts(query, fiveUTRlength, threeUTRlength)
  CDSranges = CDSranges %>% as.data.frame() %>%
    dplyr::mutate(type = 'CDS', 
                  gene_id = mcols(queryGRanges)$gene_id[1], 
                  transcript_id = mcols(queryGRanges)$transcript_id[1]) %>%
    dplyr::mutate(phase = cumsum(width%%3)%%3)
  CDSranges$phase = c(0, head(CDSranges$phase, - 1))
  CDSranges = makeGRangesFromDataFrame(CDSranges, keep.extra.columns = TRUE)
  
  return(list(ORF_considered = CDSranges))
}