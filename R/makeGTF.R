makeGTF <- function(inputGRanges, report_df, makeGTF){
  infoLog('Creating GTF file')
  if (makeGTF == TRUE){
    out.dir = getwd()
  } else{
    out.dir = makeGTF
  }
  
  # collect CDS information from reported dataframe
  report_CDS = report_df %>% 
    dplyr::select(dplyr::starts_with('ORF_considered')) %>%
    dplyr::filter(!is.na(ORF_considered.type)) %>%
    dplyr::rename_all(function(x){return(substr(x, start = 16, stop = nchar(x)))}) %>%
    dplyr::mutate(source = 'NMDer') 
  
  # collect transcript information from input GRanges
  input_transcripts = report_df %>% 
    dplyr::select(NMDer_ID, Transcript_ID) %>%
    dplyr::left_join(as.data.frame(inputGRanges), by = c('Transcript_ID' = 'transcript_id')) %>%
    dplyr::mutate(source = 'NMDer', transcript_id = NMDer_ID) %>%
    dplyr::select(seqnames:type,phase:gene_id,transcript_id)
  
  # combine transcript annotation and CDS information
  output_gtf = suppressWarnings(dplyr::bind_rows(input_transcripts, report_CDS))
  rtracklayer::export(output_gtf, paste0(out.dir, '/NMDer.gtf'), format = 'gtf')
}