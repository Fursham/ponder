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
#' 
#' 
#' 
#' 

prepareInputs <- function(queryfile, ref, fasta = NULL, in_format, ref_format) {
  
  # import assembled transcripts
  infoLog('Importing assembled transcripts...', logf, quiet)
  
  if (!file.exists(queryfile)){
    stopLog('Input transcript file do not exist')
  }
  
  # use file extension format if provided by user
  if (!is.null(in_format)) {
    fileformat = in_format
  }
  
  inputGRanges = rtracklayer::import(queryfile)
  
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
    
    # use file extension format if provided by user and import file
    if (!is.null(ref_format)) {
      fileformat = ref_format
    } 
    basicGRanges = rtracklayer::import(ref)
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
preTesting <- function(inputGRanges, basicGRanges, genome, correct_chrom, correct_gene_id, primary_gene_id, secondary_gene_id, clusters) {
  
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
                                secondary_gene_id=secondary_gene_id, 
                                clusters = clusters)
  } 
  # if user opts out of this service, program will return warning if there are unmatched gene_ids
  #   program will also update inputGRanges to include match_level and old_gene_id metadata
  else {
    
    # get a list of unmatched gene_ids
    nonstand_ids = elementMetadata(inputGRanges) %>% as.data.frame() %>%
      dplyr::select(gene_id = gene_id) %>%
      dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
      dplyr::distinct() %>%
      dplyr::mutate(status = TRUE)
    
    if (nrow(nonstand_ids) > 0) {
      
      # return warning
      warnLog(sprintf('%s transcripts have unmatched gene_ids and will not be analyzed', nrow(nonstand_ids)))
      
      # and add match_level metadata to inputGRanges
      unique_ids = elementMetadata(inputGRanges) %>% as.data.frame() %>%
        dplyr::select(gene_id = gene_id) %>% 
        dplyr::mutate(match_level = 0) %>%
        dplyr::left_join(nonstand_ids, by=c('gene_id'='gene_id')) %>%
        dplyr::mutate_at(vars(match_level), funs(ifelse(!is.na(status), 5, .))) %>%
        dplyr::select(gene_id, match_level)
      mcols(inputGRanges)$match_level = unique_ids$match_level
      mcols(inputGRanges)$old_gene_id = mcols(inputGRanges)$gene_id
      
      rm(unique_ids)
    }
    rm(nonstand_ids)
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
matchGeneIDs <- function(inputGRanges, basicGRanges, primary_gene_id=NULL, secondary_gene_id=NULL, clusters = 4) {
  
  # testing and matching gene_ids
  infoLog('Checking and matching gene_ids...')
  
  # preserve gene_ids from input GTF
  elementMetadata(inputGRanges)$old_gene_id <- elementMetadata(inputGRanges)$gene_id
  
  # get a dataframe of all gene_ids
  #  depending on the input arguments, this dataframe is constructed differently
  #  dataframe will contain headers gene_id, new_id, match_level and might contain secondary
  #     gene_id is a reference to the orignal gene_id
  #     new_id is the newly changed id
  #     match_level describes the level at which the id has been changed
  #     0 -> ids are found in ref; 
  #     1 -> ids are matched by secondary_gene_id
  #     2 -> ids are matched by appending ENS... suffix
  #     3 -> ids are matched by secondary_gene_id followed by appending ENS... suffix
  #     4 -> ids are matched by matching overlapping coordinates
  #     5 -> id could not be matched and will be skipped from analysis
  if (!is.null(primary_gene_id) & !is.null(secondary_gene_id)) {
    unique_ids = elementMetadata(inputGRanges) %>% as.data.frame() %>%
      dplyr::select(gene_id = primary_gene_id, secondary = secondary_gene_id, name = gene_name)
  } else if (!is.null(primary_gene_id)) {
    unique_ids = elementMetadata(inputGRanges) %>% as.data.frame() %>%
      dplyr::select(gene_id = primary_gene_id, name = gene_name)
  } else {
    unique_ids = elementMetadata(inputGRanges) %>% as.data.frame() %>%
      dplyr::select(gene_id = gene_id, name = gene_name)
  }
  unique_ids = unique_ids %>% 
    dplyr::mutate(new_id = gene_id, match_level = 0)
  
  # count number of non standard ID before correction
  nonstand_before = nrow(elementMetadata(inputGRanges) %>% as.data.frame() %>%
                           dplyr::select(gene_id = gene_id) %>%
                           dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
                           dplyr::distinct())
  
  # There are 3 correction steps below.
  # correction 1 and 2 will initiate only if arguments are provided
  # correction 3 will always be activated, but may give false positive output
  # note that in each step, ONLY the non-standard gene_ids will be corrected
  
  # correction 1: replace primary_gene_id with secondary_gene_id, IF both args are provided
  if (!is.null(primary_gene_id) & !is.null(secondary_gene_id)) {
    infoLog(sprintf('-> Attempting to replace %s with %s...', primary_gene_id, secondary_gene_id))
    
    # count number of non-standard ids before matching
    countsbefore = nrow(unique_ids %>%
                          dplyr::select(gene_id = new_id) %>%
                          dplyr::distinct() %>%
                          dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))
    
    # replace nonstandard primary_gene_ids from inputGRanges with its secondary_gene_ids
    unique_ids = unique_ids %>% 
      dplyr::mutate_at(vars(secondary), funs(ifelse(new_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE, ., NA))) %>%
      dplyr::mutate_at(vars(match_level), funs(ifelse(!is.na(secondary), 1, .))) %>%
      dplyr::mutate_at(vars(new_id), funs(ifelse(!is.na(secondary), as.character(secondary), .))) %>%
      dplyr::select(gene_id, new_id, name, match_level)
    
    # count number of non-standard ids after matching
    countsafter = nrow(unique_ids %>%
                         dplyr::select(gene_id = new_id) %>%
                         dplyr::distinct() %>%
                         dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))
    
    # report number of IDs corrected
    infoLog(sprintf('-> %s gene IDs matched', (countsbefore - countsafter)))
  }
  
  # correction 2: replace primary_gene_id with basic gene ID IF:
  # at least primary_gene_id is provided and if it starts with 'ENS'
  if (!is.null(primary_gene_id)) {
    infoLog('--> Attempting to match ensembl gene_ids...')
    
    
    # count number of non-standard ids after matching
    countsbefore = nrow(unique_ids %>%
                          dplyr::select(gene_id = new_id) %>%
                          dplyr::distinct() %>%
                          dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))
    
    
    # make dataframe of non_standard primary_gene_ids starting with "ENS"
    gene_id_df = unique_ids %>%
      dplyr::select(gene_id = new_id) %>%
      dplyr::distinct() %>%
      dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
      dplyr::filter(startsWith(gene_id, "ENS") == TRUE)
    
    # proceed with correction if there is an "ENS" gene_id annotation
    if (nrow(gene_id_df) > 0) {
      # make dataframe of annotated gene_id starting with "ENS"
      basic_gene_id_df = elementMetadata(basicGRanges) %>% as.data.frame() %>%
        dplyr::select(ens_id = gene_id, ens_name = gene_name) %>%
        dplyr::distinct() %>%
        dplyr::filter(startsWith(ens_id, "ENS") == TRUE)
      
      # append suffix of ENS gene_ids on both dataframes 
      # eg. ENSMUST00000029780.11 to ENSMUST00000029780
      gene_id_df$appended_id <-  sapply(strsplit(gene_id_df$gene_id, "\\."), `[`, 1)
      basic_gene_id_df$appended_id <-  sapply(strsplit(basic_gene_id_df$ens_id, "\\."), `[`, 1)
      
      # combine both data frames based on common 'appended_id' and filter lines with empty 'ens_id'
      gene_id_df = dplyr::left_join(gene_id_df, basic_gene_id_df, by = 'appended_id') %>%
        dplyr::filter(!is.na(ens_id))
      
      # replace new_id on unique_ids dataframe with ens_id, if match is found
      unique_ids = unique_ids %>% left_join(gene_id_df, by=c('new_id'='gene_id')) %>%
        dplyr::mutate_at(vars(match_level), funs(ifelse(!is.na(ens_id), as.integer(.)+1, .))) %>%
        dplyr::mutate_at(vars(new_id), funs(ifelse(!is.na(ens_id), as.character(ens_id), .))) %>%
        dplyr::mutate_at(vars(name), funs(ifelse(!is.na(ens_id), as.character(ens_name), .))) %>%
        dplyr::select(gene_id, new_id, name, match_level)
      
      # count number of non-standard ids after matching
      countsafter = nrow(unique_ids %>%
                           dplyr::select(gene_id = new_id) %>%
                           dplyr::distinct() %>%
                           dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))
      
      # report number of IDs corrected
      infoLog(sprintf('--> %s gene IDs matched', (countsbefore - countsafter)))
      
    } else {
      # check if there is any ENS ids to begin with
      ENSids = unique_ids %>%
        dplyr::distinct() %>%
        dplyr::filter(startsWith(new_id, "ENS") == TRUE)
      
      if (nrow(ENSids) > 0) {
        infoLog('--> All ensembl gene ids have been matched')
        
      } else {
        warnLog('--> No ensembl gene ids found in query')
      }
    }
  }
  
  # correction 3: correct gene_ids by finding overlapping regions.
  infoLog('---> Attempting to correct gene_ids by finding overlapping coordinates...')
  
  
  # count number of unmatched ids before matching
  countsbefore = nrow(unique_ids %>%
                        dplyr::select(gene_id = new_id) %>%
                        dplyr::distinct() %>%
                        dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))
  
  # make dataframe of non_standard gene_ids
  gene_id_df = unique_ids %>%
    dplyr::select(gene_id = new_id, test_id = gene_id) %>%
    dplyr::distinct() %>%
    dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
    as.data.frame()
  
  # this is the main worker of this script
  #   function will retrieve first gene_id from gencode that overlaps with first exon of query
  #   the below parameters give a true positive rate of 81% 
  if (nrow(gene_id_df) > 0) {
    
    group <- rep(1:clusters, length.out = nrow(gene_id_df))
    gene_id_df <- bind_cols(tibble(group), gene_id_df)
    cluster <- create_cluster(cores = clusters, quiet = TRUE)
    
    parallel_df = gene_id_df %>% partition(group, cluster = cluster)
    parallel_df %>%
      # Assign libraries
      cluster_library("GenomicRanges") %>%
      cluster_library("dplyr") %>%
      # Assign values (use this to load functions or data to each core)
      cluster_assign_value("inputGRanges", inputGRanges) %>%
      cluster_assign_value("basicGRanges", basicGRanges)
    
    gene_id_df <- parallel_df %>% # Use by_group party_df
      do(dplyr::mutate(., ref_gene_id = sapply(.$test_id, function(x) {
        
        getfirstOverlap = findOverlaps(
          inputGRanges[mcols(inputGRanges)$gene_id == x &
                         mcols(inputGRanges)$type == "exon"][1],
          basicGRanges, 
          select = "first")
        
        return(ifelse(!is.na(getfirstOverlap), 
                      mcols(basicGRanges[getfirstOverlap])$gene_id,
                      NA))
      }))) %>%
      collect() %>% # Special collect() function to recombine partitions
      as.data.frame() %>%
      dplyr::filter(!is.na(ref_gene_id))
    
    # get gene names of matched gene_ids
    gene_names_df = elementMetadata(basicGRanges) %>% as.data.frame() %>%
      dplyr::select(id = gene_id, new_name = gene_name) %>%
      dplyr::distinct()
    gene_id_df = gene_id_df %>% dplyr::left_join(gene_names_df, by = c('ref_gene_id'='id'))
    
    # replace gene_id on assembled transcript with predicted gene_id from overlap
    unique_ids = unique_ids %>% left_join(gene_id_df, by=c('new_id'='gene_id')) %>%
      dplyr::mutate_at(vars(match_level), funs(ifelse(!is.na(ref_gene_id), 4, .))) %>%
      dplyr::mutate_at(vars(new_id), funs(ifelse(!is.na(ref_gene_id), as.character(ref_gene_id), .))) %>%
      dplyr::mutate_at(vars(name), funs(ifelse(!is.na(ref_gene_id), as.character(new_name), .))) %>%
      dplyr::select(gene_id, new_id, name, match_level)
    
    # count number of unmatched ids after matching
    countsafter = nrow(unique_ids %>%
                         dplyr::select(gene_id = new_id) %>%
                         dplyr::distinct() %>%
                         dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE))
    
    # report number of IDs corrected
    infoLog(sprintf('---> %s gene IDs matched', (countsbefore - countsafter)))
    
  } else {
    infoLog('---> Skipped, all gene IDs matched')
  }
  
  
  # count remaining number of non_standard gene_ids and number of corrected ids
  nonstand_id_list = unique_ids %>%
    dplyr::select(gene_id = new_id) %>%
    dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
    dplyr::distinct() %>%
    dplyr::mutate(status = TRUE)
  nonstand_after = nrow(nonstand_id_list)
  corrected_ids = nonstand_before - nonstand_after
  
  if (nonstand_after > 0) {
    
    # if there are remaining unmatched gene_ids, update match_level to 5.
    #   these transcripts will be skipped from further analysis
    unique_ids = unique_ids %>% left_join(nonstand_id_list, by=c('new_id'='gene_id')) %>%
      dplyr::mutate_at(vars(match_level), funs(ifelse(!is.na(status), 5, .))) %>%
      dplyr::select(gene_id, new_id, match_level)
  }
  
  # update gene_id and add match_level to inputGRanges
  mcols(inputGRanges)$gene_id = unique_ids$new_id
  mcols(inputGRanges)$gene_name = unique_ids$name
  mcols(inputGRanges)$match_level = unique_ids$match_level
  
  # report pre-testing analysis and return inputGRanges
  infoLog(sprintf('Total gene_ids corrected: %s', corrected_ids))
  infoLog(sprintf('Remaining number of non-standard gene_ids: %s', nonstand_after))
  if (nonstand_after > 0) {
    warnLog('Transcripts with non-standard gene_ids will be skipped from analysis')
    
  } 
  
  # cleanup and return
  rm(list = c('unique_ids','nonstand_id_list'))
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
                  Chrom = seqnames,
                  Strand = strand) %>%
    dplyr::mutate(NMDer_ID = as.character(NA),
                  Tx_coordinates = as.character(NA),
                  ORF_considered = as.character(NA),
                  Alt_tx = as.logical(NA),
                  annotatedStart = as.logical(NA),
                  predictedStart = as.logical(NA))
  
  # prepare df of gencode_basic transcripts
  basicTX_df = elementMetadata(basicGRanges) %>% as.data.frame() %>%
    dplyr::filter(type == 'transcript') %>%
    dplyr::select(Gene_ID = gene_id, Ref_TX_ID = transcript_id) %>%
    dplyr::distinct() %>% as.data.frame()
  
  # join report_df and basicTX_df and generate a full list of query transcripts to be compared with every reference transcript
  # and add in NMDerIDs
  combined_report_df = dplyr::right_join(basicTX_df, report_df, by = 'Gene_ID') %>%
    dplyr::select(NMDer_ID, Gene_ID, Original_Gene_ID, Match_level, 
                  Gene_Name, Transcript_ID, Ref_TX_ID, Chrom, Strand, 
                  Tx_coordinates, annotatedStart, predictedStart, Alt_tx,
                  ORF_considered) %>% 
    dplyr::mutate(rownum = rownames(.)) %>%
    rowwise() %>% dplyr::mutate_at(vars(NMDer_ID), 
                                   funs(paste("NMDer", formatC(as.integer(rownum), width=7, flag="0"), sep=""))) %>%
    dplyr::select(-rownum)
  
  # prepare databases
  inputDB = makeTxDbFromGRanges(inputGRanges)
  basicDB = makeTxDbFromGRanges(basicGRanges)
  
  # prepare exon structures by transcripts/CDSs
  inputExonsbyTx = exonsBy(inputDB, by="tx", use.names=TRUE)
  basicExonsbyCDS = cdsBy(basicDB, by="tx", use.names=TRUE)
  basicExonsbyTx = exonsBy(basicDB, by="tx", use.names=TRUE)
  
  # cleanup unused variables
  rm(list = c('inputDB','basicDB'))
  return(list(combined_report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx))
}





