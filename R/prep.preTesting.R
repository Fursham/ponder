#' NMDer workflow: Check and match input annotations
#'
#' @description Check the chromosome ID and gene ID of query, 
#' reference and genome objects and if requested, attempts to match these metadata 
#'
#' @param inputGRanges Query GRanges object
#' @param basicGRanges Reference GRanges object
#' @param genome Genome sequence as a Biostring object
#' @param correct_chrom If TRUE, attempt to match chromosome ID
#' @param correct_gene_id If TRUE, attempt to match gene ID. This is crudely done by intersecting query transcript coordinates with reference. To increase chance of true matching, user may provide 
#' @param primary_gene_id See 'correct_gene_id' for details
#' @param secondary_gene_id See 'correct_gene_id' for details
#'
#' @return Query GRanges object that have been matched

preTesting <- function(inputGRanges, basicGRanges, genome, correct_chrom, correct_gene_id, primary_gene_id, secondary_gene_id) {
  
  # testing and correcting chromosome names on query and annotated transcripts
  infoLog('Checking and correcting chromosome names...', logf, quiet)
  
  if (correct_chrom == TRUE & length(slotNames(genome)) != 5) {
    if (seqlevelsStyle(inputGRanges) != seqlevelsStyle(genome)) {
      newStyle <- mapSeqlevels(seqlevels(inputGRanges), seqlevelsStyle(genome))
      newStyle = newStyle[!is.na(newStyle)]
      inputGRanges <- renameSeqlevels(inputGRanges, newStyle)
      
      if (any(!seqlevels(inputGRanges)%in%seqlevels(genome))) {
        seqlevels(inputGRanges, pruning.mode = 'tidy') <- as.vector(newStyle)
        warnLog('Non-standard chromosome IDs in query were removed')
      }
    }
    if (seqlevelsStyle(basicGRanges) != seqlevelsStyle(genome)) {
      newStyle <- mapSeqlevels(seqlevels(basicGRanges), seqlevelsStyle(genome))
      newStyle = newStyle[!is.na(newStyle)]
      basicGRanges <- renameSeqlevels(basicGRanges, newStyle)
      
      if (any(!seqlevels(basicGRanges)%in%seqlevels(genome))) {
        seqlevels(basicGRanges, pruning.mode = 'tidy') <- as.vector(newStyle)
        warnLog('Non-standard chromosome IDs in reference were removed')
      }
    }
  } 
  # if user opts out of this correction service, program will test if there are non-standard IDs and return a warning
  # program will continue
  else {
    if (length(slotNames(genome)) == 5) {
      warnLog('Please ensure that user-provided fasta file contain common seqlevels as query and reference')
    }
    
    if (any(!seqlevels(inputGRanges)%in%seqlevels(genome))) {
      warnLog('Non-standard chromosome IDs in query were found')
    }
    if (any(!seqlevels(basicGRanges)%in%seqlevels(genome))) {
      warnLog('Non-standard chromosome IDs in reference were found')
    }
  }
  
  
  # correct gene ids if requested by user
  #   else, calculate the number of unmatched gene_ids and return warning
  if (correct_gene_id == TRUE) {
    inputGRanges = matchGeneIDs(inputGRanges, basicGRanges, 
                                primary_gene_id=primary_gene_id, 
                                secondary_gene_id=secondary_gene_id)
  } 
  # if user opts out of this service, program will return warning if there are unmatched gene_ids
  #   program will also update inputGRanges to include match_level and old_gene_id metadata
  else {
    
    # prepare a df with a list of gene_ids found in reference
    basicGRanges.genelist = basicGRanges %>% as.data.frame() %>%
      dplyr::select(gene_id) %>%
      dplyr::distinct() %>%
      dplyr::mutate(matched = TRUE)   # this column functions to annotate matched genes later on using join
    
    
    # get a list of unmatched gene_ids
    inputGRanges = suppressMessages(inputGRanges %>% as.data.frame() %>%
                                      dplyr::left_join(basicGRanges.genelist) %>%
                                      dplyr::mutate(match_level = ifelse(is.na(matched), 5, 0)) %>%
                                      dplyr::mutate(old_gene_id = gene_id) %>%
                                      dplyr::select(-matched))
    
    count_nonstand_ids = inputGRanges %>% 
      dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
      dplyr::filter(match_level == 5) %>%
      nrow()
    
    if (count_nonstand_ids > 0) {
      
      # return warning
      warnLog(sprintf('%s transcripts have unmatched gene_ids and will not be analyzed', count_nonstand_ids))
    }
    
    # convert inputGRanges back to GRanges object
    inputGRanges = makeGRangesFromDataFrame(inputGRanges, keep.extra.columns = TRUE)
  }
  return(list('inputGRanges' = inputGRanges, 
              'basicGRanges' = basicGRanges))
}
