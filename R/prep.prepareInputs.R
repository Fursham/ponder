#' NMDer workflow: Import files
#' 
#' @description 
#' Imports mandatory files for NMDer analysis
#' 
#' @param queryFile Path to query GTF/GFF3 transcript annotation file
#' @param refFile Path to ref GTF/GFF3 transcript annotation file
#' @param fasta BSGenome object (preferred) or path to fasta file
#' @param user_query_format optional argument to specify the query annotation format ('gtf','gff3')
#' @param user_ref_format optional argument to specify the refFileerence annotation format ('gtf','gff3')
#' 
#' @return
#' A list containing (1) GRanges object for query, (2) GRanges object for 
#' ref annotation, and (3) DNAstring containing genome sequence
#' 

prepareInputs <- function(queryFile, refFile, fasta, 
                          user_query_format, user_ref_format) {
  
  # import assembled transcripts
  infoLog('Importing query transcript file', logf, quiet)
  
  if (!file.exists(queryFile)){
    stopLog('Query transcript file do not exist')
  }
  
  # check file extension 
  query_format = tools::file_ext(queryFile)

  # return if infile is a txt file with no user_query_format input
  if (query_format == 'txt' & is.null(user_query_format)) {
    stopLog('Query file contain .txt extension but query_format argument not provided')
  } else if (query_format == 'txt' & !is.null(user_query_format)) {
    query_format = user_query_format
  }

  # process a gff3 query file to be compatible for analysis
  if (query_format == 'gff3') {
    inputGRanges = rtracklayer::import(queryFile, format = 'gff3')
    
    # convert gff3 to gtf
    inputGRanges = gff3togtfconvert(inputGRanges, comment = 'query')
    
  } else if (query_format == 'gtf') {
    inputGRanges = rtracklayer::import(queryFile, format = 'gtf')
  } else {
    stopLog('Query file format not supported')
  }

  if ('*' %in% (BiocGenerics::strand(inputGRanges) %>% levels())) {
    warnLog('Query file contain transcripts with no strand information. These will be removed')
    inputGRanges = inputGRanges[BiocGenerics::strand(inputGRanges) != '*']
  }
  
  
  # load provided genome_basic assembly, or import user refFileerence annotation
  if (is(refFile, 'GenomicRanges')) {
    infoLog(sprintf('Loading gencode_basic transcripts'), logf, quiet)
    
    basicGRanges = refFile
  } else if (is.character(refFile)){
    infoLog('Importing user-provided reference transcript annotations', logf, quiet)
    
    # return if file does not exist
    if (!file.exists(refFile)){
      stopLog('Reference annotation file do not exist')
    }
    
    # check file extension 
    ref_format = tools::file_ext(refFile)
    
    # return if infile is a txt file with no user_query_format input
    if (ref_format == 'txt' & is.null(user_ref_format)) {
      stopLog('Reference file contain .txt extension but reference_format argument not provided')
    } else if (ref_format == 'txt' & !is.null(user_ref_format)) {
      ref_format = user_ref_format
    }
    
    
    # process a gff3 query file to be compatible for analysis
    if (ref_format == 'gff3') {
      basicGRanges = rtracklayer::import(refFile, format = 'gff3')
      
      if (!'gene' %in% basicGRanges$type) {
        stopLog('Please ensure that reference GFF3 file contain gene entries')
      }
      basicGRanges = gff3togtfconvert(basicGRanges, comment = 'reference')
      
    } else if (ref_format == 'gtf') {
      basicGRanges = rtracklayer::import(refFile, format = 'gtf')
    } else {
      stopLog('Reference file format not supported')
    }
    
  }
  if ('*' %in% (BiocGenerics::strand(basicGRanges) %>% levels())) {
    warnLog('Reference annotation contain transcripts with no strand information. These will be removed')
    basicGRanges = basicGRanges[BiocGenerics::strand(basicGRanges) != '*']
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

