
#' Title
#'
#' @param df 
#' @param output_dir 
#'
#' @return
#' @export
#'
#' @examples
generateGTF <- function(df, output_dir) {
  
  infoLog('Writing GTF file...', logf, quiet)
  
  outfile = sprintf("%s/NMDer.gtf", output_dir)
  write ('#', file = outfile)
  gtfdf = apply(df, 1, function(y) {
    
    tempdf = c()
    
    # get line info and extract exon coords
    exon_coords = unlist(strsplit(as.character(y[['Tx_coordinates']]), ';')) %>% as.data.frame() %>% 
      tidyr::separate('.', into = c('start', 'end'), sep = '-') %>%
      dplyr::mutate(exon_num = row_number())
    
    # prepare transcript info and add to tempdf
    txstart = ifelse(y[['Strand']]== '-', exon_coords[nrow(exon_coords),1], exon_coords[1,1])
    txstop = ifelse(y[['Strand']]== '-',exon_coords[1,2], exon_coords[nrow(exon_coords),2])
    
    
    txmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s; is_nmd %s;', 
                         dQuote(y[['Gene_ID']]), dQuote(y[['NMDer_ID']]), dQuote(y[['Gene_Name']]), y[['is_NMD']])
    tx = list(chrom = as.character(y[['Chrom']]), source = 'NMDer', type = 'transcript', start = as.integer(txstart), end = as.integer(txstop),
              score = 1000, strand = as.character(y[['Strand']]), frame = '.', meta = txmetainfo)
    tempdf= rbind(tempdf, tx)
    
    
    # prepare exon info and add to tempdf
    
    out = apply(exon_coords, 1, function(x) {
      exmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s; exon_number %s;', 
                           dQuote(y[['Gene_ID']]), dQuote(y[['NMDer_ID']]), dQuote(y[['Gene_Name']]), x[['exon_num']])
      ex = list(chrom = as.character(y[['Chrom']]), source = 'NMDer', type = 'exon', start = as.integer(x[['start']]), end = as.integer(x[['end']]),
                score = 1000, strand = as.character(y[['Strand']]), frame = '.', meta = exmetainfo)
      return(ex)
      
      #write(paste(ex, collapse = '\t'), file = outfile, append = TRUE)
    })
    out = do.call(rbind, out)
    tempdf = rbind(tempdf, out)
    
    
    # return cds info if transcript has one
    if (!is.na(y[['ORF_considered']])) {
      cds_coords = unlist(strsplit(as.character(y[['ORF_considered']]), ';')) %>% as.data.frame() %>% 
        tidyr::separate('.', into = c('start', 'end'), sep = '-') %>%
        dplyr::mutate(exon_num = row_number())
      
      # extract information about frames of the CDS
      cdsIRanges = IRanges(start = as.integer(cds_coords[,1]), end = as.integer(cds_coords[,2]))
      frames = c(0, head(cumsum(width(cdsIRanges)%%3)%%3, -1))
      cds_coords$frame = frames
      
      # prepare and add start_codon information
      startcodonstart = ifelse(y[['Strand']]== '-', as.integer(cds_coords[1,2]) - 2, cds_coords[1,1])
      startcodonend = ifelse(y[['Strand']]== '-', cds_coords[1,2], as.integer(cds_coords[1,1])+2)
      startmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s;', 
                              dQuote(y[['Gene_ID']]), dQuote(y[['NMDer_ID']]), dQuote(y[['Gene_Name']]))
      start = list(chrom = as.character(y[['Chrom']]), source = 'NMDer', type = 'start_codon', start = startcodonstart, end = startcodonend,
                   score = 1000, strand = as.character(y[['Strand']]), frame = as.character(0), meta = startmetainfo)
      tempdf= rbind(tempdf, start)
      
      #write(paste(start, collapse = '\t'), file = outfile, append = TRUE)
      
      
      
      # prepare and add stop_codon information
      stopcodonstart = ifelse(y[['Strand']]== '-', cds_coords[nrow(cds_coords),1], as.integer(cds_coords[nrow(cds_coords),2]) - 2)
      stopcodonend = ifelse(y[['Strand']]== '-', as.integer(cds_coords[nrow(cds_coords),1])+2, cds_coords[nrow(cds_coords),2])
      
      stopframe = frames[length(frames)]
      stopmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s;',
                             dQuote(y[['Gene_ID']]), dQuote(y[['NMDer_ID']]), dQuote(y[['Gene_Name']]))
      stop = list(chrom = as.character(y[['Chrom']]), source = 'NMDer', type = 'stop_codon', start = stopcodonstart, end = stopcodonend,
                  score = 1000, strand = as.character(y[['Strand']]), frame = as.character(stopframe), meta = stopmetainfo)
      tempdf= rbind(tempdf, stop)
      #write(paste(stop, collapse = '\t'), file = outfile, append = TRUE)
      
      
      # prepare CDS cooords information
      out = apply(cds_coords, 1, function(x) {
        cdsmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s;', 
                              dQuote(y[['Gene_ID']]), dQuote(y[['NMDer_ID']]), dQuote(y[['Gene_Name']]))
        cds = list(chrom = as.character(y[['Chrom']]), source = 'NMDer', type = 'CDS', start = as.integer(x[['start']]), end = as.integer(x[['end']]),
                   score = 1000, strand = as.character(y[['Strand']]), frame = as.character(x[['frame']]), meta = cdsmetainfo)
        return(cds)
        
        #write(paste(cds, collapse = '\t'), file = outfile, append = TRUE)
      })
      out = do.call(rbind, out)
      tempdf = rbind(tempdf, out)
      
    }
    
    return(tempdf)
  })
  gtfdf = do.call(rbind, gtfdf)
  write.table(gtfdf, file = outfile, row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE, append = TRUE)
}



filterdata <- function(df) {
  
  # simplify dataframe by taking NMDer_ID, Transcript_ID and coverages values and each transcript by decreasing coverage
  newdf = dplyr::select(df, NMDer_ID, Transcript_ID, Shared_coverage) %>% as.data.frame() %>%
    dplyr::arrange(Transcript_ID, desc(Shared_coverage)) %>%
    dplyr::distinct(Transcript_ID, .keep_all = TRUE)
  
  # filter df
  outputdf = dplyr::filter(df, NMDer_ID%in%newdf$NMDer_ID == TRUE)
  return(outputdf)
}

summSplicing <- function(df) {
  
  # extract columns from dataframe corresponding to splicing classes
  splicing_df = df %>% as.data.frame() %>% 
    dplyr::select(NMDer_ID, Gene_ID, Gene_Name, Transcript_ID, Ref_TX_ID,
                  Chrom, Strand, is_NMD, CE:NMDcausing) %>% 
    dplyr::mutate_at(vars(CE:ir), funs(ifelse(is.na(.), 0, lengths((strsplit(as.character(.), ';'))))))
  
  return(splicing_df)
}
