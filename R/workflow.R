#' PONDER workflow: Prepare PONDER analysis
#'
#' @description Import transcript annotation file, 
#' match chromosome levels and gene IDs and prepare for NMD prediction analysis
#'
#' @param query Mandatory. Path to query GTF/GFF3 transcript annotation file
#' @param reference Mandatory. Path to reference GTF/GFF3 transcript annotation file. 
#' @param fasta Mandatory. BSGenome object (preferred) or path to fasta file
#' 
#' @param query_format Optional argument to specify the query annotation format ('gtf','gff3'). 
#' Mandatory if query contains '.txt' extension filename
#' @param reference_format Optional argument to specify the reference annotation format ('gtf','gff3'). 
#' Mandatory if reference contains '.txt' extension filename
#' 
#' @param match_chrom Supplementary feature. If TRUE, program will attempt to match chromosome names of 
#' query and reference to fasta genome to ensure consistent naming across input files.
#' @param match_geneIDs Supplementary feature to attempt to match gene IDs in query file
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
#' 
#' @param primary_gene_id See match_geneIDs argument
#' @param secondary_gene_id See match_geneIDs argument
#'
#' @return S4 object containing dataframes and objects for downstream NMD prediction analysis
#' @export
#'
#' @examples 
#' library("BSgenome.Mmusculus.UCSC.mm10")
#' preppedObject = prepNMDer(testQuery, testRef, Mmusculus, match_geneIDs = TRUE)
#' preppedObject = prepNMDer(testQuery, testRef, Mmusculus, match_geneIDs = TRUE, primary_gene_id = 'gene_id', secondary_gene_id = 'ref_gene_id')
#' 
prePonder <- function(query,
                 reference, 
                 fasta,
                 query_format = NULL,
                 reference_format = NULL,
                 match_chrom = FALSE,
                 match_geneIDs = FALSE,
                 primary_gene_id = NULL,
                 secondary_gene_id = NULL) {

  options(warn=-1)
  # check for mandataory arguments
  if (any(missing(query), missing(reference), missing(fasta))) {
    stopLog('Missing mandatory arguments')
  } 
  
  # import and/or load query file(s)
  unpack[inputGRanges, basicGRanges, genome] = 
    prepareInputs(query, reference, fasta, query_format, reference_format)
  
  # matching chromosome names and gene IDs
  unpack[inputGRanges, basicGRanges] = preTesting(inputGRanges, basicGRanges, genome, 
                            match_chrom, match_geneIDs, 
                            primary_gene_id, secondary_gene_id)
  
  # prepare dataframes and transcript GRanges objects
  unpack[report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx] = 
    prepareAnalysis(inputGRanges, basicGRanges)
  
  options(warn=0)
  
  # pack dataframes and biostring objects into an S4 object for return
  preppedOutput = packNMDer(df = report_df, 
                            inputGRanges = inputGRanges, 
                            inputTranscripts = inputExonsbyTx, 
                            basicTranscripts = basicExonsbyTx, 
                            basicCDS = basicExonsbyCDS, 
                            fasta = genome)
  
  return(preppedOutput)
}



#' PONDER workflow: Run PONDER analysis
#'
#' @description Execute core analysis
#'
#' @param prepObject 
#' S4 object from prepNMDer output
#' 
#' @param testNMD
#' If TRUE, program will test transcripts for primary NMD feature
#'  
#' @param testOtherFeatures 
#' If TRUE, program will test transcripts for other NMD features
#' 
#' @param PTC_dist 
#' Minimum distance of PTC to last exon-junction to trigger NMD.
#' Default: 50bp
#' 
#' @param testAS 
#' If TRUE, program will test for alternatively spliced segments
#' 
#' @param makeGTF 
#' If TRUE, program will output a GTF file of query
#' Provide path to output
#' 
#' @param clusters 
#' Number of cores to run the program. Default = 4
#'
#' @return df with analysis
#' @export
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#'
#' @examples
#' 
runPonder <- function(prepObject,
                 testNMD = TRUE,
                 testOtherFeatures = FALSE,
                 PTC_dist = 50,
                 testAS = FALSE,
                 makeGTF = FALSE,
                 clusters = 4) {
  
  options(warn=-1)
  
  # unpack custom S4 object
  unpack[report_df, inputGRanges, inputExonsbyTx, basicExonsbyTx, basicExonsbyCDS, genome] = 
    unPack(prepObject)

  
  # separate
  report_df = report_df %>%
    dplyr::mutate(Coverage = 0,
                  ORF_considered = as.character(NA),
                  ORF_start = as.character('Not found'),
                  ORF_found = FALSE)
  
  report_df_unmatched = report_df %>%
    dplyr::filter(Gene_match_level == 5) %>%
    dplyr::mutate(Ref_transcript_ID = as.character(Ref_transcript_ID)) %>%
    dplyr::select(-ORF_considered)
  
  report_df = report_df %>% 
    dplyr::filter(Gene_match_level != 5)
  
  infoLog('Preparing clusters for analysis')
  cluster = parallel::makeCluster(clusters)
  doParallel::registerDoParallel(cluster)
  
  infoLog('Predicting NMD features')
  report_df <- foreach::foreach(i=1:nrow(report_df), .combine = dplyr::bind_rows) %dopar% {
  runMain(report_df[i,], inputExonsbyTx, basicExonsbyCDS,
                          basicExonsbyTx, genome, testNMD, 
                       PTC_dist,testOtherFeatures, testAS)
  } %>%
    as.data.frame() %>%
    dplyr::arrange(NMDer_ID) %>% dplyr::select(-ORF_considered)
  
  doParallel::stopImplicitCluster()
  
  # OLD parallel code using multidplyr 
  #
  #
  # # run analysis using single core or multidplyr
  # if(clusters == 1){
  #   infoLog('Predicting NMD features')
  #   report_df <- report_df %>% 
  #     dplyr::do(runMain(., inputExonsbyTx, basicExonsbyCDS,
  #                basicExonsbyTx, genome, testNMD, 
  #                PTC_dist,testOtherFeatures, testAS)) %>%
  #     as.data.frame() %>%
  #     dplyr::arrange(NMDer_ID) %>% dplyr::select(-ORF_considered)
  # } else {
  #   # run NMD analysis in parallel
  #   infoLog('Preparing clusters for analysis')
  #   group <- rep(1:clusters, length.out = nrow(report_df))
  #   report_df <- dplyr::bind_cols(tibble::tibble(group), report_df)
  #   cluster <- multidplyr::new_cluster(clusters)
  #   
  #   parallel_df = report_df %>% 
  #     dplyr::group_by(group) %>% 
  #     multidplyr::partition(cluster)
  #   cluster %>%
  #     # Assign libraries
  #     multidplyr::cluster_library("NMDer") %>%
  #     multidplyr::cluster_library("Biostrings") %>%
  #     multidplyr::cluster_library("BSgenome") %>%
  #     multidplyr::cluster_library("dplyr") %>%
  #     multidplyr::cluster_assign("basicExonsbyCDS" = basicExonsbyCDS) %>%
  #     multidplyr::cluster_assign("inputExonsbyTx" = inputExonsbyTx) %>%
  #     multidplyr::cluster_assign("basicExonsbyTx" = basicExonsbyTx) %>% 
  #     multidplyr::cluster_assign("genome" = genome) %>%
  #     multidplyr::cluster_assign("testNMD" = testNMD) %>%
  #     multidplyr::cluster_assign("PTC_dist" = PTC_dist) %>%
  #     multidplyr::cluster_assign("testOtherFeatures" = testOtherFeatures) %>%
  #     multidplyr::cluster_assign("testAS" = testAS)
  #   
  #   infoLog('Predicting NMD features')
  #   report_df <- parallel_df %>% # Use by_group party_df
  #     do(runMain(., inputExonsbyTx, basicExonsbyCDS,
  #                basicExonsbyTx, genome, testNMD, 
  #                PTC_dist,testOtherFeatures, testAS)) %>%
  #     collect() %>% # Special collect() function to recombine partitions
  #     as.data.frame() %>%
  #     dplyr::arrange(NMDer_ID) %>% dplyr::select(-group, -ORF_considered)
  # }

  # combine report with unmatched entries
  output_df = report_df %>%
    dplyr::select(-dplyr::starts_with('ORF_considered')) %>%
    dplyr::bind_rows(report_df_unmatched) %>%
    dplyr::distinct(NMDer_ID, .keep_all = TRUE) %>%
    dplyr::arrange(NMDer_ID)
  
  # prepare GTF transcript anntotation output
  if (makeGTF != FALSE){
    makeGTF(inputGRanges, report_df, makeGTF)
  }
  
  infoLog('Done!')
  options(warn=0)

  return(output_df)
}



