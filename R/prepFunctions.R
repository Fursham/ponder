#' _ workflow: Load and/or import input files
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW. 
#' This function will import and/or load input data for the workflow and
#' will return a list containing GRanges objects for assembled transcript and 
#' reference annotation, and fasta file for genome sequence
#' 
#' @param file filename or /path/to/file of transcript assemblies. File can be of
#' gtf, gff3 or bed format, and its filename should contain the corresponding extension
#' @param gencode reference annotation of coding transcripts for each gene family. This program
#' contain in-built gencode_basic annotation for GRCh38_hg38 and GRCm38_mm10 assemblies. These
#' assemblies can be loaded by parsing "hg38" or "mm10" into this argument. If users wish to 
#' import other assemblies, please provide filename or /path/to/file of reference CDS assemblies.
#' File can be of gtf, gff3 or bed format, and its filename should contain the corresponding extension
#' @param fasta genome sequence in fasta format. This program contain in-built fasta files for 
#' GRCh38_hg38 and GRCm38_mm10 genome. These sequences can be loaded by parsing "hg38" or "mm10" 
#' into this argument. If users wish to import other genomes, please provide filename or 
#' /path/to/file of fasta files.
#'
#' @return
#' A list containing (1) GRanges object for assembled transcripts, (2) GRanges object for 
#' reference CDS assembly, and (3) DNAstring containing genome sequence
#' @export
#' @import rtracklayer
#' @import Biostrings
#' @import tools
#' 
#' 
#' 
#' 

prepareInputs <- function(queryfile, ref, fasta = NULL, user_query_format, user_ref_format) {
  
  # import assembled transcripts
  infoLog('Importing assembled transcripts...', logf, quiet)
  
  if (!file.exists(queryfile)){
    stopLog('Input transcript file do not exist')
  }
  
  # check file extension 
  query_format = tools::file_ext(queryfile)

  # return if infile is a txt file with no user_query_format input
  if (query_format == 'txt' & is.null(user_query_format)) {
    stopLog('Input transcript file contain .txt extension. Please provide transcript annotation format using argument query_format')
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
    stopLog('Input file format not supported')
  }
  
  if ('*'%in%strand(inputGRanges)) {
    warnLog('Query GTF file contain transcripts with no strand information. These will be removed')
    inputGRanges = inputGRanges[strand(inputGRanges) != '*']
  }
  

  
  # load provided genome_basic assembly, or import user reference annotation
  if (is(ref, 'GenomicRanges')) {
    infoLog(sprintf('Loading gencode_basic transcripts...'), logf, quiet)
    
    basicGRanges = ref
  } else if (is.character(ref)){
    infoLog('Importing user-provided reference annotations...', logf, quiet)
    
    # return if file does not exist
    if (!file.exists(ref)){
      stopLog('Reference annotation file do not exist')
    }
    
    
    
    
    basicGRanges = rtracklayer::import(ref)
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


#' _ workflow: Testing and correcting transcript assembly
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW. 
#' This function will test and correct the input transcript assembly for consistent chromosome names 
#' and gene ids with CDS assembly and genome sequence.
#' 
#' Non-consistent Gene IDs are corrected by replacing it with Gene IDs from reference CDS assembly
#' that contain overlapping exons. However, this strategy may be inaccurate and may return 
#' false-postive gene ids. To improve the correction, we provide two upper layers of correction
#' that can be initiated by providing the mandatory 'primary_gene_id' and 'secondary_gene_id'
#' arguments
#' 
#' Primary_gene_id refer to the name of the metadata column containing gene ids for each assemnbled transcript. 
#' This is typically 'gene_id'
#' Secondary_gene_id refer to the name of the metadata column that may contain a reference gene id
#' For example: 
#' NEED TO DESCRIBE
#'  
#' But if both args are provided, 
#' function will try to replace 
#' 
#' How to check for metadata columns
#'
#' @param inputGRanges 
#' GRanges object containing transcript and exon features. It is recommended to create this object
#' by importing gtf/gff3 files using rtracklayer::import
#' @param basicGRanges 
#' #' GRanges object containing reference transcript,exon and CDS features. 
#' It is recommended to create this object by importing gtf/gff3 files using rtracklayer::import
#' @param genome
#' DNAstring object containing genome sequence
#' @param primary_gene_id 
#' Optional argument. Typically, this is 'gene_id'. Refer to description for details
#' @param secondary_gene_id 
#' Optional argument. Refer to description for details
#'
#' @return
#' Function will return a modified inputGRanges object
#' 
#' @export
#' @import GenomeInfoDb
#'
#' @examples
preTesting <- function(inputGRanges, basicGRanges, genome, correct_chrom, correct_gene_id, primary_gene_id, secondary_gene_id) {
  
  # testing and correcting chromosome names on query and annotated transcripts
  infoLog('Checking and correcting chromosome names...', logf, quiet)
  
  if (correct_chrom == TRUE) {
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
      dplyr::filter(type == 'transcript', match_level == 5) %>%
      nrow()

    if (count_nonstand_ids > 0) {
      
      # return warning
      warnLog(sprintf('%s transcripts have unmatched gene_ids and will not be analyzed', count_nonstand_ids))
    }
    
    # convert inputGRanges back to GRanges object
    inputGRanges = makeGRangesFromDataFrame(inputGRanges, keep.extra.columns = TRUE)
  }
  return(inputGRanges)
}




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
#' Assembled transcripts in query. Argument can be provided as name of GTF/GFF3 file in current working directory,
#' or path_to_file. Argument can also be provided as a GenomicRanges object containing transcript and/or exon ranges.
#' @param ref 
#' Annotated reference transcripts. This program contain gencode basic annotations from mm10 and hg38 assemblies  
#' that will be loaded using 'gencode_basics[['mm10']] or gencode_basics[['hg38']] as input to argument. Argument can also be provided as name of 
#' GTF/GFF3 file in current working directory, or path_to_file.
#' @param query_format 
#' Format of query file (GTF/GFF3). If this argument is not provided and input to query argument
#' is a filename, program will attempt to read file extension.
#' @param ref_format 
#' Format of reference file (GTF/GFF3). If this argument is not provided and input to query argument
#' is a filename, program will attempt to read file extension.
#' @param primary_gene_id 
#' Name of the primary gene id in query file. Input to this argument is typically 'gene_id'
#' @param secondary_gene_id 
#' Name of the secondary gene id in query file. Example of input to this arguement is 'ref_gene_id'
#' 
#' @param outputfile 
#' Name of output file (Default: 'matched_geneIDs.gtf'). File will be saved in current working directory
#'
#' @return
#' GTF file with matched Gene IDs if makefile == TRUE. 
#' GRanges object if makefile == FALSE.
#' 
#' @export
#'
#' @examples
#' 
#' matchGeneIDs(testData, gencode_basics[['mm10']], primary_gene_id = 'gene_id', secondary_gene_id = 'ref_gene_id', clusters = 4)
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
  inputGRanges = suppressMessages(inputGRanges %>% as.data.frame() %>%
    dplyr::mutate(old_gene_id = gene_id, match_level = 0) %>% 
    dplyr::left_join(basicGRanges.genelist))
  
  
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
                                        group_by(transcript_id) %>%
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
                                        ungroup())
      
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
    #   and match those gene_names to the reference gene_names
    #   and change the match_level of matched transcripts
    unmatched_df = inputGRanges %>%
      dplyr::filter(is.na(matched), exon_number == 1) %>%
      dplyr::select(seqnames, start, end, strand, transcript_id)
    unmatched_granges = makeGRangesFromDataFrame(unmatched_df, keep.extra.columns = TRUE)
    
    matched_df = mergeByOverlaps(unmatched_granges, basicGRanges) %>% 
      as.data.frame() %>%
      dplyr::select(transcript_id, basic_gene_id = gene_id, basic_gene_name = gene_name) %>%
      dplyr::distinct(transcript_id, .keep_all = TRUE)
    
    inputGRanges = suppressMessages(inputGRanges %>%
                                      dplyr::select(-matched) %>%
                                      dplyr::left_join(matched_df) %>%
                                      dplyr::mutate(gene_id = ifelse(!is.na(basic_gene_id),
                                                                     basic_gene_id, 
                                                                     gene_id)) %>%
                                      dplyr::mutate(gene_name = ifelse(!is.na(basic_gene_name),
                                                                       basic_gene_name, 
                                                                       gene_name)) %>%
                                      dplyr::mutate(match_level = ifelse(!is.na(basic_gene_id),
                                                                         4, 
                                                                         match_level)) %>%
                                      dplyr::select(-basic_gene_id, -basic_gene_name) %>%
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
  
  inputGRanges = makeGRangesFromDataFrame(inputGRanges, keep.extra.columns = TRUE)
  return(inputGRanges)
}


#' _ workflow: Preparing for analysis
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW.
#' This function will prepare output dataframes and GRanges objects for the _Program 
#'
#' @param inputGRanges 
#' GRanges object containing transcript and exon features. It is recommended to create this object
#' by importing gtf/gff3 files using rtracklayer::import
#' @param basicGRanges 
#' #' GRanges object containing reference transcript,exon and CDS features. 
#' It is recommended to create this object by importing gtf/gff3 files using rtracklayer::import
#'
#' @return
#' Function will return a list containing:
#' (1) dataframe of a list of assembled transcripts and its information for output
#' (2) dataframe of a list of reference CDS transcripts
#' (3) GRangesList of exon coordinates of assembled transcripts grouped by transcript name
#' (4) GRangesList of exon coordinates of reference CDS grouped by transcript name
#' 
#' @export
#' @import dplyr
#' @import GenomicFeatures
#'
#' @examples
prepareAnalysis <- function(inputGRanges, basicGRanges) {
  
  infoLog('Preparing databases, transcripts and output files...', logf, quiet)
  
  # prepare output dataframe by getting transcript information from inputGRanges
  report_df = inputGRanges %>% as.data.frame() %>%
    dplyr::filter(type == 'transcript') %>%
    dplyr::select(Gene_ID = gene_id, 
                  Original_Gene_ID = old_gene_id,
                  Gene_Name = gene_name,
                  Match_level = match_level,
                  Transcript_ID = transcript_id,
                  seqnames,
                  start,
                  end,
                  Strand = strand) %>%
    dplyr::mutate(NMDer_ID = as.character(NA),
                  Alt_tx = as.logical(NA),
                  annotatedStart = as.logical(NA),
                  predictedStart = as.logical(NA)) %>%
    tidyr::unite(Tx_coord., c(start, end), sep = '-') %>%
    tidyr::unite(Tx_coord, c(seqnames, Tx_coord.), sep = ':') 
  
  # prepare df of gencode_basic transcripts
  basicTX_df = basicGRanges %>% as.data.frame() %>%
    dplyr::filter(type == 'transcript') %>%
    dplyr::select(Gene_ID = gene_id, Ref_TX_ID = transcript_id) %>%
    dplyr::distinct()
  
  # join report_df and basicTX_df and generate a full list of query transcripts to be compared with every reference transcript
  # and add in NMDerIDs
  combined_report_df = suppressMessages(report_df %>%
      dplyr::left_join(basicTX_df) %>%
      dplyr::select(NMDer_ID, Gene_ID, Original_Gene_ID, Match_level, 
                  Gene_Name, Transcript_ID, Ref_TX_ID, Tx_coord, 
                  Strand, annotatedStart, predictedStart, Alt_tx)) 

  # prepare databases
  inputDB = makeTxDbFromGRanges(inputGRanges)
  basicDB = makeTxDbFromGRanges(basicGRanges)
  
  # prepare exon structures by transcripts/CDSs
  inputExonsbyTx = exonsBy(inputDB, by="tx", use.names=TRUE) %>% as.data.frame()
  basicExonsbyCDS = cdsBy(basicDB, by="tx", use.names=TRUE) %>% as.data.frame()
  basicExonsbyTx = exonsBy(basicDB, by="tx", use.names=TRUE) %>% as.data.frame()
  
  # remove comparisons in which the reference transcript do not have a CDS
  basicCDS = basicExonsbyCDS %>% as.data.frame() %>%
    dplyr::select(Ref_TX_ID = group_name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(CDS = TRUE)
  
  combined_report_df = suppressMessages(combined_report_df %>%
    dplyr::left_join(basicCDS) %>%
    dplyr::filter((!is.na(CDS) & !is.na(Ref_TX_ID)) | Match_level == 5) %>%
    dplyr::select(-CDS))
  
  # combine query transcripts with multiple comparisons to reference
  combined_report_df = combined_report_df %>% 
    group_by(Transcript_ID) %>%
    mutate(Ref_TX_ID = list(as.character(Ref_TX_ID))) %>%
    ungroup() %>%
    distinct(.keep_all = TRUE) %>%
    dplyr::mutate(NMDer_ID = paste0('NMDer', formatC(as.integer(row_number()), width=7, flag='0')))

  
  # cleanup unused variables
  rm(list = c('inputDB','basicDB'))
  return(list(combined_report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx))
}





