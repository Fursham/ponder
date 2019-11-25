getCDSranges <- function(query, fiveUTRlength, threeUTRlength){
  
  CDSranges = resizeTranscripts(query, fiveUTRlength, threeUTRlength)
  CDSranges = CDSranges %>% as.data.frame() %>%
    dplyr::mutate(type = 'CDS', 
                  gene_id = mcols(query)$gene_id[1], 
                  transcript_id = mcols(query)$transcript_id[1]) %>%
    dplyr::mutate(phase = cumsum(width%%3)%%3) %>%
    dplyr::select(seqnames:end, strand, type, phase, gene_id, transcript_id)
  CDSranges$phase = c(0, head(CDSranges$phase, - 1))
  CDSranges = GenomicRanges::makeGRangesFromDataFrame(CDSranges, keep.extra.columns = TRUE)
  
  return(list(ORF_considered = CDSranges))
}