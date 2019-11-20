#' NMDer workflow: Import files
#' 
#' @description 
#' Imports mandatory files for NMDer analysis
#' 
#' @param queryfile Name of the query GTF/GFF3 transcript annotation file
#' @param ref Name of the reference GTF/GFF3 transcript annotation file.
#' Alternatively, user may choose to use mm10 or hg38 gencode basic annotation that 
#' comes pre-loaded with NMDer
#' @param fasta Genome sequence in the form of Biostrings object (preferred) 
#' or name of fasta genome sequence file for import
#' @param user_query_format optional argument to specify the query annotation format ('gtf','gff3')
#' @param user_ref_format optional argument to specify the reference annotation format ('gtf','gff3')
#' 
#' @return
#' A list containing (1) GRanges object for query, (2) GRanges object for 
#' reference annotation, and (3) DNAstring containing genome sequence
#' 
#' @import rtracklayer
#' @import Biostrings
#' @import tools
#' 
#' 


prepareInputs <- function(queryfile, ref, fasta = NULL, 
                          user_query_format, user_ref_format) {
  
  # import assembled transcripts
  infoLog('Importing query transcript file', logf, quiet)
  
  if (!file.exists(queryfile)){
    stopLog('Query transcript file do not exist')
  }
  
  # check file extension 
  query_format = tools::file_ext(queryfile)
  
  # return if infile is a txt file with no user_query_format input
  if (query_format == 'txt' & is.null(user_query_format)) {
    stopLog('Query file contain .txt extension but query_format argument not provided')
  } else if (query_format == 'txt' & !is.null(user_query_format)) {
    query_format = user_query_format
  }
  
  # process a gff3 query file to be compatible for analysis
  # basically, gff3 annotation contain different headers as compared to gtf2
  # so this function attempts to extract the 1) transcript_id, 2) gene_id, 
  # 3)gene_name and 4)exon_number of each transcript, which are important information
  # for downstream functions
  if (query_format == 'gff3') {
    inputGRanges = rtracklayer::import(queryfile, format = 'gff3')
    
    # removes line of type 'gene', extract transcript_id from ID (for mRNA/transcript types)
    # or from Parent header (the other types), add gene_id and gene_name to every entry
    inputGRanges = inputGRanges %>% as.data.frame() %>%
      dplyr::filter(type != 'gene') %>%
      dplyr::mutate(Parent = as.character(Parent)) %>%
      dplyr::mutate(transcript_id = ifelse(type %in% c('mRNA', 'transcript'), ID, Parent))%>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(desc(width)) %>% 
      dplyr::mutate(gene_id = Parent[1]) %>%
      dplyr::mutate(gene_name = Name[1]) %>%
      ungroup()
    
    # duplicate dataframe and filter out 'exon' type entries, sort and add exon_number
    inputGRanges.exons = inputGRanges %>%
      dplyr::filter(type == 'exon') %>%
      dplyr::group_by(transcript_id) %>% 
      dplyr::arrange(ifelse(strand == '-', desc(start), start)) %>%
      dplyr::mutate(exon_number = row_number()) %>% 
      ungroup()
    
    # combine both dataframes and sort by transcript_id
    inputGRanges = suppressMessages(inputGRanges %>%
                                      left_join(inputGRanges.exons) %>% 
                                      arrange(transcript_id))
    
    inputGRanges = makeGRangesFromDataFrame(inputGRanges, keep.extra.columns = TRUE)
    
  } else if (query_format == 'gtf') {
    inputGRanges = rtracklayer::import(queryfile, format = 'gtf')
  } else {
    stopLog('Query file format not supported')
  }
  
  if ('*'%in%strand(inputGRanges)) {
    warnLog('Query file contain transcripts with no strand information. These will be removed')
    inputGRanges = inputGRanges[strand(inputGRanges) != '*']
  }
  
  
  
  # load provided genome_basic assembly, or import user reference annotation
  if (is(ref, 'GenomicRanges')) {
    infoLog(sprintf('Loading gencode_basic transcripts'), logf, quiet)
    
    basicGRanges = ref
  } else if (is.character(ref)){
    infoLog('Importing user-provided reference transcript annotations', logf, quiet)
    
    # return if file does not exist
    if (!file.exists(ref)){
      stopLog('Reference annotation file do not exist')
    }
    
    # check file extension 
    ref_format = tools::file_ext(ref)
    
    # return if infile is a txt file with no user_query_format input
    if (ref_format == 'txt' & is.null(user_ref_format)) {
      stopLog('Reference file contain .txt extension but reference_format argument not provided')
    } else if (ref_format == 'txt' & !is.null(user_ref_format)) {
      ref_format = user_ref_format
    }
    
    
    # process a gff3 query file to be compatible for analysis
    # basically, gff3 annotation contain different headers as compared to gtf2
    # so this function attempts to extract the 1) transcript_id, 2) gene_id, 
    # 3)gene_name and 4)exon_number of each transcript, which are important information
    # for downstream functions
    if (ref_format == 'gff3') {
      basicGRanges = rtracklayer::import(ref, format = 'gff3')
      
      # removes line of type 'gene', extract transcript_id from ID (for mRNA/transcript types)
      # or from Parent header (the other types), add gene_id and gene_name to every entry
      basicGRanges = basicGRanges %>% as.data.frame() %>%
        dplyr::filter(type != 'gene') %>%
        dplyr::mutate(Parent = as.character(Parent)) %>%
        dplyr::mutate(transcript_id = ifelse(type %in% c('mRNA', 'transcript'), ID, Parent))%>%
        dplyr::group_by(transcript_id) %>%
        dplyr::arrange(desc(width)) %>% 
        dplyr::mutate(gene_id = Parent[1]) %>%
        dplyr::mutate(gene_name = Name[1]) %>%
        ungroup()
      
      # duplicate dataframe and filter out 'exon' type entries, sort and add exon_number
      basicGRanges.exons = basicGRanges %>%
        dplyr::filter(type == 'exon') %>%
        dplyr::group_by(transcript_id) %>% 
        dplyr::arrange(ifelse(strand == '-', desc(start), start)) %>%
        dplyr::mutate(exon_number = row_number()) %>% 
        ungroup()
      
      # combine both dataframes and sort by transcript_id
      basicGRanges = suppressMessages(basicGRanges %>%
                                        left_join(basicGRanges.exons) %>% 
                                        arrange(transcript_id))
      
      basicGRanges = makeGRangesFromDataFrame(basicGRanges, keep.extra.columns = TRUE)
      
    } else if (ref_format == 'gtf') {
      basicGRanges = rtracklayer::import(ref, format = 'gtf')
    } else {
      stopLog('Reference file format not supported')
    }
    
  }
  if ('*'%in%strand(basicGRanges)) {
    warnLog('Reference annotation contain transcripts with no strand information. These will be removed')
    basicGRanges = basicGRanges[strand(basicGRanges) != '*']
  }
  
  
  # import genome sequence
  if (is(fasta, 'BSgenome') | typeof(fasta) == 'S4') {
    infoLog('Loading genome sequence', logf, quiet)
    genome = fasta 
  } 
  else if (is.character(fasta)) {
    infoLog('Importing fasta file...', logf, quiet)
    
    if (!file.exists(fasta)){
      stopLog('Fasta file do not exist')
    } else {
      genome = Biostrings::readDNAStringSet(fasta, format="fasta")
    }
  } else {
    genome = NULL
  }
  return(list('inputGRanges' = inputGRanges, 
              'basicGRanges' = basicGRanges, 
              'genome' = genome))
}

