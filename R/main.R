#' Title
#'
#' @param query 
#' @param reference 
#' @param fasta 
#' @param query_format 
#' @param reference_format 
#' @param match_chrom 
#' @param match_geneIDs 
#' @param primary_gene_id 
#' @param secondary_gene_id 
#' @param clusters 
#'
#' @return
#' @export
#'
#' @examples
prepNMDer <- function(query,
                 reference, 
                 fasta,
                 query_format = NULL,
                 reference_format = NULL,
                 match_chrom = FALSE,
                 match_geneIDs = FALSE,
                 primary_gene_id = NULL,
                 secondary_gene_id = NULL,
                 clusters = 4) {

  options(warn=-1)
  # check for mandataory arguments
  if (any(missing(query), missing(reference), missing(fasta))) {
    stopLog('Missing mandatory arguments')
  } 
  
  # import and/or load query file(s)
  unpack[inputGRanges, basicGRanges, genome] = 
    prepareInputs(query, reference, fasta,query_format, reference_format)
  
  # unpacking objects
  #inputGRanges = packedInput$inputGRanges
  #basicGRanges = packedInput$basicGRanges
  #genome = packedInput$genome
  
  # matching chromosome names and gene IDs
  inputGRanges = preTesting(inputGRanges, basicGRanges, genome, 
                            match_chrom, match_geneIDs, 
                            primary_gene_id, secondary_gene_id,
                            clusters)
  
  # prepare dataframes and transcript GRanges objects
  unpack[report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx] = 
    prepareAnalysis(inputGRanges, basicGRanges, output_dir)
  
  options(warn=0)
  return(list(report_df, inputExonsbyTx, basicExonsbyTx, basicExonsbyCDS, genome))
}


#' Title
#'
#' @param prepObject 
#' @param testNMD 
#' @param other_features 
#' @param PTC_dist 
#' @param testAS 
#' @param clusters 
#'
#' @return
#' @import multidplyr
#' @export
#'
#' @examples
runNMDer <- function(prepObject,
                 testNMD = TRUE,
                 testOtherFeatures = FALSE,
                 PTC_dist = 50,
                 testAS = FALSE,
                 clusters = 4) {
  
  options(warn=-1)
  
  # unpack object
  unpack[report_df, inputExonsbyTx, basicExonsbyTx, basicExonsbyCDS, genome] = 
    prepObject
  
  infoLog('Running NMDer')

  # run NMD analysis in parallel
  group <- rep(1:clusters, length.out = nrow(report_df))
  report_df <- bind_cols(tibble(group), report_df)
  cluster <- create_cluster(cores = clusters, quiet = TRUE)
  
  parallel_df = report_df %>% partition(group, cluster = cluster)
  parallel_df %>%
    # Assign libraries
    cluster_library("NMDer") %>%
    cluster_library("Biostrings") %>%
    cluster_library("BSgenome") %>%
    cluster_library("dplyr") %>%
    cluster_assign_value("basicExonsbyCDS", basicExonsbyCDS) %>%
    cluster_assign_value("inputExonsbyTx", inputExonsbyTx) %>%
    cluster_assign_value("basicExonsbyTx", basicExonsbyTx) %>% 
    cluster_assign_value("genome", genome) %>%
    cluster_assign_value("testNMD", testNMD) %>%
    cluster_assign_value("PTC_dist", PTC_dist) %>%
    cluster_assign_value("testOtherFeatures", testOtherFeatures) %>%
    cluster_assign_value("testAS", testAS)
  
  report_df <- parallel_df %>% # Use by_group party_df
    do(runMain(., inputExonsbyTx, basicExonsbyCDS,
                       basicExonsbyTx, genome, testNMD, 
                       PTC_dist,testOtherFeatures, testAS)) %>%
    collect() %>% # Special collect() function to recombine partitions
    as.data.frame() %>%
    dplyr::arrange(NMDer_ID) %>% dplyr::select(-group)
  
  infoLog('Done!')
  options(warn=0)

  return(report_df)
}



matchIDs <- function(query,
                     reference,
                     query_format = NULL,
                     reference_format = NULL,
                     primary_gene_id = NULL,
                     secondary_gene_id = NULL,
                     outputfile = 'matched_geneIDs.gtf',
                     clusters = 4) {
  
  options(warn=-1)
  # check for mandataory arguments
  if (any(missing(query), missing(reference), missing(fasta))) {
    stopLog('Missing mandatory arguments')
  } 
  
  # import and/or load query file(s)
  unpack[inputGRanges, basicGRanges] = 
    prepareInputs(query, reference, query_format, reference_format)

  # matching chromosome names and gene IDs
  inputGRanges = matchGeneIDs(inputGRanges, basicGRanges, primary_gene_id, 
                              secondary_gene_id, clusters = 4)
  
  # export
  rtracklayer::export(inputGRanges, outputfile)
  infoLog(sprintf('Done. GTF saved as %s in current working directory', outputfile))
  
  options(warn=0)
  return(list(report_df, inputExonsbyTx, basicExonsbyTx, basicExonsbyCDS, genome))
}
