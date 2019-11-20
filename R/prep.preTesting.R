#' NMDer workflow: Check and match input annotations
#'
#' @description This program checks the chromosome ID and gene ID of query, 
#' reference and genome objects and if requested, attempts to match these entries 
#'
#' @param inputGRanges Query GRanges object
#' @param basicGRanges Reference GRanges object
#' @param genome Genome sequence as a Biostring object
#' @param correct_chrom Supplementary feature. If TRUE, program will attempt to match chromosome names of 
#' query and reference to fasta genome to ensure consistent naming across input files.
#' @param correct_gene_id Supplementary feature to attempt to match gene IDs in query file
#' to reference file. This is key in grouping query transcripts to reference gene families for comparison.
#' 
#' 
#' Matching is done at three levels with increasing accuracy:
#' 
#' 1. Crudely intersecting query coordinates with reference. Invoked by setting match_geneIDs to TRUE
#' 
#' 2. Trim ensembl-style gene IDs and attempt matching. Invoked by providing name of 
#' gene ID header (typically 'gene_id') from gtf file to primary_gene_id argument
#' 
#' 3. Replace query gene ID with a secondary gene ID and attempt matching. Invoked by providing name of 
#' secondary gene ID header (for example 'ref_gene_id') from gtf file to secondary_gene_id argument 
#' @param primary_gene_id See 'correct_gene_id' for details
#' @param secondary_gene_id See 'correct_gene_id' for details
#'
#' @return Query GRanges object that have been matched

preTesting <- function(inputGRanges, basicGRanges, genome, correct_chrom, correct_gene_id, primary_gene_id, secondary_gene_id) {
  
  # testing and correcting chromosome names on query and annotated transcripts
  infoLog('Checking and correcting chromosome names...', logf, quiet)
  
  if (correct_chrom == TRUE & length(slotNames(genome)) != 5) {
    
    # attempt to match input and reference GRanges to genome
    inputGRanges = matchSeqLevels(inputGRanges, genome)
    basicGRanges = matchSeqLevels(basicGRanges, genome)
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
