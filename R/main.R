#' Run _ Program
#' 
#' @description
#' 
#'
#' @param query I
#' @param reference 
#' @param fasta 
#' @param primary_gene_id 
#' @param secondary_gene_id 
#' @param other_features 
#' @param makegtf 
#' @param PTC_dist 
#' @param output_dir
#' @param quiet
#'
#' @return
#' @export
#' @import zeallot
#' @import multidplyr
#' 
#' @examples
#' library("BSgenome.Mmusculus.UCSC.mm10")
#' run(testData, mm10, BSgenome.Mmusculus.UCSC.mm10, match_geneIDs = TRUE, primary_gene_id = 'gene_id', secondary_gene_id = 'ref_gene_id')
#' 
#' 
#' 
#' 
run <- function(query,
                reference, 
                fasta,
                query_format = NULL,
                reference_format = NULL,
                match_chrom = FALSE,
                match_geneIDs = FALSE,
                primary_gene_id = NULL,
                secondary_gene_id = NULL, 
                other_features = FALSE,
                filterbycoverage = TRUE,
                make_gtf = TRUE,
                PTC_dist = 50,
                output_dir = "NMDer",
                clusters = 4,
                quiet = FALSE) {
  
  # create output directory and logfile
  dir.create(output_dir, showWarnings = FALSE)
  assign('logf', file(sprintf('%s/NMDer.log', output_dir), 'wt'), envir = .GlobalEnv)
  assign("quiet", quiet, envir = .GlobalEnv)
  options(warn=-1)
  
  # check for mandataory arguments
  if (any(missing(query), missing(reference), missing(fasta))) {
    stopLog('Missing mandatory arguments', logf)
  } 

  # import and/or load query file(s)
  c(inputGRanges, basicGRanges, genome) %<-% prepareInputs(query, reference, 
                                                           fasta, query_format, 
                                                           reference_format)

  # test for standard chromosome names, and gene_ids between query files
  inputGRanges = preTesting(inputGRanges, basicGRanges, genome, 
                            match_chrom, match_geneIDs, 
                            primary_gene_id, secondary_gene_id,
                            clusters)
  
  # prepare dataframe and databases
  c(report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx) %<-% 
    prepareAnalysis(inputGRanges, basicGRanges, output_dir)
  
  # cleanup
  rm(list = c('reference','fasta','inputGRanges','basicGRanges'))
  
  # run NMD analysis in parallel
  infoLog('Detecting NMD features...', logf, quiet)
  
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
    cluster_assign_value("PTC_dist", PTC_dist) %>%
    cluster_assign_value("other_features", other_features) %>%
    cluster_assign_value("logf", logf) %>%
    cluster_assign_value("quiet", quiet)
  
  report_df <- parallel_df %>% # Use by_group party_df
    do(testNMDfeatures(., inputExonsbyTx, basicExonsbyCDS,
                       basicExonsbyTx,genome, PTC_dist,other_features)) %>%
    collect() %>% # Special collect() function to recombine partitions
    as.data.frame() %>%
    dplyr::arrange(NMDer_ID)
  rm(parallel_df)

    
  # prepare outputs
  outputAnalysis(report_df, filterbycoverage, other_features, make_gtf, 
                 output_dir, query, reference, PTC_dist)
  gc()
}

