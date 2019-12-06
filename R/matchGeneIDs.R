#' Match Gene IDs from query GTF/GFF3 file
#' 
#' @description 
#' This function will match and correct Gene IDs from a query assembled transcript file, using a 
#' transcript annotation as reference. 
#' 
#' The default approach to this correction relies on finding overlaps between transcripts in query with
#' transcripts in reference. Using this method alone could result in false positive matches (19 percent false positives).
#' To improve this, users have an option to activate two additional layers of matching.
#' (1) Matching by ENSEMBL Gene_IDs. If both query and reference transcript annotations containg Ensembl-style
#' Gene IDs, this program will try to match both IDs in a less stringent manner. This correction can be invoked
#' by providing the 'primary_gene_id' argument
#' 
#' (2) Matching by secondary Gene_IDs. Depending on the transcript assembly program, GTF/GFF3 annotations 
#' may contain additional comments on the transcript information. This may include a distinct 
#' secondary Gene ID annotation that potentially matches with the reference. To invoke this correction,
#' provide 'primary_gene_id' and 'secondary_gene_id' arguments. To determine if your transcript assembly contain 
#' possible secondary Gene IDs, try importing query GTF file using rtracklayer package and check its metadata 
#' columns
#' 
#' 
#' @param query 
#' Mandatory. Path to query GTF/GFF3 transcript annotation file
#' @param ref 
#' Mandatory. Path to reference GTF/GFF3 transcript annotation file. 
#' @param primary_gene_id 
#' Name of the primary gene id in query file. Input to this argument is typically 'gene_id'
#' @param secondary_gene_id 
#' Name of the secondary gene id in query file. Example of input to this arguement is 'ref_gene_id'
#' 
#' @return
#' Gene_id-matched query GRanges
#' 
#' @examples
#' 
#' 
#' 
matchGeneIDs <- function(inputGRanges, basicGRanges, primary_gene_id=NULL, secondary_gene_id=NULL) {
  
  ###################################################################################################################
  # this function will attempt to match gene_ids from input to reference
  #   there are 2 optional matching functions and 1 constitutive matching function
  #   the optional sub-function can be invoked by providing primary_gene_id and secondary_gene_id (matching 1)
  #   or by providing primary_gene_id only (matching 2)
  #   depending on the the degree of matching done, a match_level is assigned and the number represents:
  #     0 -> ids are found in ref; 
  #     1 -> ids are matched by secondary_gene_id
  #     2 -> ids are matched by appending ENS... suffix
  #     3 -> ids are matched by secondary_gene_id followed by appending ENS... suffix
  #     4 -> ids are matched by matching overlapping coordinates
  #     5 -> id could not be matched and will be skipped from analysis
  ###################################################################################################################
  
  # testing and matching gene_ids
  infoLog('Checking and matching gene_ids...')
  
  # prepare a df with a list of gene_ids found in reference
  basicGRanges.genelist = basicGRanges %>% as.data.frame() %>%
    dplyr::select(gene_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(matched = TRUE)   # this column functions to annotate matched genes later on using join
  
  # convert input GRanges object to dataframe for parsing and adding meta information
  inputGRanges = inputGRanges %>% as.data.frame() %>%
    dplyr::mutate(old_gene_id = gene_id, match_level = 0) %>% 
    dplyr::left_join(basicGRanges.genelist, by = 'gene_id')
  
  
  # count number of non standard ID before correction
  nonstand_before = inputGRanges %>% 
    dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
    dplyr::filter(is.na(matched)) %>%
    nrow()
  
  # Matching function 1: replace primary_gene_id with secondary_gene_id, IF both args are provided
  if (!is.null(primary_gene_id) & !is.null(secondary_gene_id)) {
    
    # count number of non-standard ids before matching
    countsbefore = inputGRanges %>% 
      dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
      dplyr::filter(is.na(matched)) %>%
      nrow()
    
    
    # proceed with matching only if there are unmatched gene ids
    if (countsbefore > 0){
      
      infoLog(sprintf('-> Attempting to correct gene ids by replacing %s with %s...', primary_gene_id, secondary_gene_id))
      
      # prepare a df with a list of gene_ids found in reference
      basicGRanges.genelist.1 = basicGRanges %>% as.data.frame() %>%
        dplyr::select(gene_id) %>%
        dplyr::distinct() %>%
        dplyr::mutate(matched = TRUE)
      
      # core of the function. 
      #   this will replace the primary_gene_id with the secondary_gene_id (if present)
      #   and change the match_level of matched transcripts to 1
      inputGRanges = suppressMessages(inputGRanges %>% 
                                        dplyr::mutate(!!primary_gene_id := ifelse(is.na(matched) & !is.na(get(secondary_gene_id)), 
                                                                                  get(secondary_gene_id), 
                                                                                  get(primary_gene_id))) %>%
                                        dplyr::mutate(match_level = ifelse(is.na(matched) & !is.na(get(secondary_gene_id)), 
                                                                           1, 
                                                                           match_level)) %>%
                                        dplyr::select(-matched) %>%
                                        dplyr::left_join(basicGRanges.genelist.1))
      
      # count number of non-standard ids after matching
      countsafter = inputGRanges %>% 
        dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
        dplyr::filter(is.na(matched)) %>%
        nrow()
      
      # report number of IDs corrected
      infoLog(sprintf('-> %s transcripts corrected for gene ids', (countsbefore - countsafter)))
      
    }
    
    
    
    
  }
  
  # Matching function 2: replace primary_gene_id with basic gene ID IF:
  # at least primary_gene_id is provided and if it starts with 'ENS'
  if (!is.null(primary_gene_id)) {

    # count number of non-standard ids before matching
    countsbefore = inputGRanges %>% 
      dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
      dplyr::filter(is.na(matched)) %>%
      nrow()
    
    # proceed with matching only if there are unmatched gene ids
    if (countsbefore > 0){
      
      infoLog('--> Attempting to match ensembl gene_ids...')
      # prepare a df with a list of reference gene_ids and append ENSEMBL style ids to remove suffixes
      basicGRanges.genelist.2 = basicGRanges %>% as.data.frame() %>%
        dplyr::select(gene_id) %>%
        dplyr::distinct() %>%
        dplyr::rowwise() %>%
        dplyr::mutate(matched = TRUE) %>%
        dplyr::mutate(appended_ens_id = ifelse(startsWith(gene_id, 'ENS'), strsplit(gene_id, split = '\\.')[[1]][1], NA)) %>%
        dplyr::filter(!is.na(appended_ens_id)) %>%
        dplyr::select(appended_ens_id, basic_gene_id = gene_id, matched)
      
      # core of the function. 
      #   this function will append ENSEMBL style primary_gene_ids to remove suffixes
      #   and match those gene_ids to the appended reference gene_ids
      #   and change the match_level of matched transcripts
      inputGRanges = suppressMessages(inputGRanges %>%
                                        dplyr::group_by(transcript_id) %>%
                                        dplyr::mutate(appended_ens_id = ifelse(startsWith(get(primary_gene_id), 'ENS') & is.na(matched), 
                                                                               strsplit(get(primary_gene_id), split = '\\.')[[1]][1], 
                                                                               as.character(NA))) %>%
                                        dplyr::select(-matched) %>%
                                        dplyr::left_join(basicGRanges.genelist.2) %>%
                                        dplyr::mutate(!!primary_gene_id := ifelse(!is.na(basic_gene_id),
                                                                                  basic_gene_id, 
                                                                                  get(primary_gene_id))) %>%
                                        dplyr::mutate(match_level = ifelse(!is.na(basic_gene_id),
                                                                           match_level + 2, 
                                                                           match_level)) %>%
                                        dplyr::select(-appended_ens_id, -basic_gene_id) %>%
                                        dplyr::ungroup())
      
      # count number of non-standard ids after matching
      countsafter= inputGRanges %>% 
        dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
        dplyr::filter(is.na(matched)) %>%
        nrow()
      
      # print out statistics of the match
      #   or print out warning if none of the genes were matched
      if (countsbefore > countsafter) {
        infoLog(sprintf('--> %s transcripts corrected for gene ids', (countsbefore - countsafter)))
      } else {
        anyEnsid = inputGRanges %>% 
          dplyr::select(gene_id) %>%
          dplyr::distinct() %>%
          dplyr::filter(startsWith(gene_id, 'ENS')) %>%
          nrow() > 0
        
        if (anyEnsid == TRUE) {
          infoLog('--> All ensembl gene ids have been matched')
        } else{
          warnLog('--> No ensembl gene ids found in query')
        }
      }
      
    }
  }
  
  # Matching function 3: correct gene_ids by finding overlapping regions.
  
  
  # count number of unmatched ids before matching
  countsbefore = inputGRanges %>% 
    dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
    dplyr::filter(is.na(matched)) %>%
    nrow()
  
  if (countsbefore == 0) {
    infoLog('--> All gene ids have been matched')
  } else {
    
    infoLog('---> Attempting to match gene_ids by finding overlapping coordinates...')
    # core of the function. 
    #   this function will gather all unmatched transcripts
    #   and attempt to find its overlap with the reference
    #   and match those gene_ids to the reference gene_ids
    #   and change the match_level of matched transcripts
    unmatched_df = inputGRanges %>%
      dplyr::filter(is.na(matched), exon_number == 1) %>%
      dplyr::select(seqnames, start, end, strand, transcript_id)
    unmatched_granges = GenomicRanges::makeGRangesFromDataFrame(unmatched_df, keep.extra.columns = TRUE)
    
    matched_df = IRanges::mergeByOverlaps(unmatched_granges, basicGRanges) %>% 
      as.data.frame() %>%
      dplyr::select(transcript_id, basic_gene_id = gene_id) %>%
      dplyr::distinct(transcript_id, .keep_all = TRUE)
    
    inputGRanges = suppressMessages(inputGRanges %>%
                                      dplyr::select(-matched) %>%
                                      dplyr::left_join(matched_df) %>%
                                      dplyr::mutate(gene_id = ifelse(!is.na(basic_gene_id),
                                                                     basic_gene_id, 
                                                                     gene_id)) %>%
                                      dplyr::mutate(match_level = ifelse(!is.na(basic_gene_id),
                                                                         4, 
                                                                         match_level)) %>%
                                      dplyr::select(-basic_gene_id) %>%
                                      dplyr::left_join(basicGRanges.genelist))
    
    # count number of unmatched ids after matching
    countsafter = inputGRanges %>% 
      dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
      dplyr::filter(is.na(matched)) %>%
      nrow()
    
    
    # report statistics of the match
    infoLog(sprintf('---> %s transcripts corrected for gene ids', (countsbefore - countsafter)))
    
    
  }
  
  
  ## post matching function
  # annotate the match_level on unmatched gene_ids
  # and cleanup the dataframe
  
  inputGRanges = inputGRanges %>%
    dplyr::mutate(match_level = ifelse(is.na(matched),
                                       5,match_level)) %>%
    dplyr::select(-matched)
    
  if('gene_name' %in% names(inputGRanges) & 'gene_name' %in% names(S4Vectors::mcols(basicGRanges))){
    basicGRanges.genelist.1 = basicGRanges %>% as.data.frame() %>%
      dplyr::select(gene_id, ref_gene_name = gene_name) %>%
      dplyr::distinct() 
    
    inputGRanges = inputGRanges %>%
      dplyr::left_join(basicGRanges.genelist.1, by = 'gene_id') %>%
      dplyr::mutate(gene_name = ifelse(match_level != 5 & is.na(gene_name),
                                         ref_gene_name,gene_name)) %>%
      dplyr::select(-ref_gene_name)
  } else {
    inputGRanges = inputGRanges %>%
      dplyr::mutate(gene_name = NA)
  } 
  
  # report pre-testing analysis and return inputGRanges
  nonstand_after = inputGRanges %>% 
    dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
    dplyr::filter(match_level == 5) %>%
    nrow()
  corrected_ids = nonstand_before - nonstand_after
  
  infoLog(sprintf('Total gene_ids corrected: %s', corrected_ids))
  infoLog(sprintf('Remaining number of non-standard gene_ids: %s', nonstand_after))
  if (nonstand_after > 0) {
    warnLog('Transcripts with non-standard gene_ids will be skipped from analysis')
    
  } 
  
  inputGRanges = GenomicRanges::makeGRangesFromDataFrame(inputGRanges, keep.extra.columns = TRUE)
  return(inputGRanges)
}


