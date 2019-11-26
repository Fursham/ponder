getORFranges <- function(query, fiveUTRlength, threeUTRlength){
  
  # prepare output list
  output = list('ORF_considered' = NA)
  
  # resize query GRanges to ORF and renew metadata info
  CDSranges = resizeGRangesTranscripts(query, fiveUTRlength, threeUTRlength)
  CDSranges = CDSranges %>% as.data.frame() %>%
    dplyr::mutate(type = 'CDS', 
                  gene_id = mcols(query)$gene_id[1], 
                  transcript_id = mcols(query)$transcript_id[1]) %>%
    dplyr::mutate(phase = cumsum(width%%3)%%3) %>%
    dplyr::select(seqnames:end, strand, type, phase, gene_id, transcript_id)
  CDSranges$phase = c(0, head(CDSranges$phase, - 1))
  output$ORF_considered  = GenomicRanges::makeGRangesFromDataFrame(CDSranges, keep.extra.columns = TRUE)

  return(output)
}