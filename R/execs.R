#' NMDer workflow: Prepare NMDer analysis
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
#' Matching is done at three levels with increasing accuracy:
#' 1. Crudely intersecting query coordinates with reference. Invoked by setting match_geneIDs to TRUE
#' 2. Trim ensembl-style gene IDs and attempt matching. Invoked by providing name of 
#' gene ID header (typically 'gene_id') from gtf file to primary_gene_id argument
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
#' preppedObject = prepNMDer(testData, mm10, Mmusculus, match_geneIDs = TRUE)
#' preppedObject = prepNMDer(testData, mm10, Mmusculus, match_geneIDs = TRUE, primary_gene_id = 'gene_id', secondary_gene_id = 'ref_gene_id')
#' 
prepNMDer <- function(query,
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
                 makeGTF = FALSE,
                 clusters = 4) {
  
  options(warn=-1)
  
  # unpack custom S4 object
  unpack[report_df, inputGRanges, inputExonsbyTx, basicExonsbyTx, basicExonsbyCDS, genome] = 
    unPack(prepObject)

  
  # separate
  report_df_unmatched = report_df %>%
    dplyr::filter(Match_level == 5) %>%
    dplyr::mutate(Ref_TX_ID = as.character(Ref_TX_ID))
  
  report_df = report_df %>% 
    dplyr::filter(Match_level != 5)

  # run NMD analysis in parallel
  infoLog('Preparing clusters for analysis')
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
  
  infoLog('Predicting NMD features')
  report_df <- parallel_df %>% # Use by_group party_df
    do(runMain(., inputExonsbyTx, basicExonsbyCDS,
                       basicExonsbyTx, genome, testNMD, 
                       PTC_dist,testOtherFeatures, testAS)) %>%
    collect() %>% # Special collect() function to recombine partitions
    as.data.frame() %>%
    dplyr::arrange(NMDer_ID) %>% dplyr::select(-group, -ORF_considered)
  
  
  output_df = report_df %>% 
    dplyr::select(-starts_with('ORF_considered')) %>%
    dplyr::bind_rows(report_df_unmatched) %>%
    dplyr::distinct(NMDer_ID, .keep_all = TRUE) %>%
    dplyr::arrange(NMDer_ID)
  
  if (makeGTF != FALSE){
    infoLog('Creating GTF file')
    if (makeGTF == TRUE){
      out.dir = getwd()
    } else{
      out.dir = makeGTF
    }
    
    report_CDS = report_df %>% 
      dplyr::select(starts_with('ORF_considered')) %>%
      dplyr::filter(!is.na(ORF_considered.type))
    names(report_CDS) = substr(names(report_CDS), start = 16, stop = nchar(report_CDS))
    report_CDS$source = 'NMDer'

    input_transcripts = report_df %>% 
      dplyr::select(NMDer_ID, Transcript_ID) %>%
      dplyr::left_join(as.data.frame(inputGRanges), by = c('Transcript_ID' = 'transcript_id')) %>%
      dplyr::mutate(source = 'NMDer', transcript_id = NMDer_ID) %>%
      dplyr::select(seqnames:type,phase:gene_id,transcript_id)
    
    output_gtf = bind_rows(input_transcripts, report_CDS)
    rtracklayer::export(output_gtf, paste0(out.dir, '/NMDer.gtf'), format = 'gtf')
  }
  
  
  
  infoLog('Done!')
  options(warn=0)

  return(output_df)
}
