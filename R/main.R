#' Run _ Program
#' 
#' @description
#' 
#'
#' @param input I
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
#' 
#' @examples
#' library("BSgenome.Mmusculus.UCSC.mm10")
#' run(testData, 'mm10', BSgenome.Mmusculus.UCSC.mm10, match_geneIDs = TRUE, primary_gene_id = 'gene_id', secondary_gene_id = 'ref_gene_id')
#' 
#' 
#' 
#' 
run <- function(input,
                reference, 
                fasta,
                input_format = NULL,
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
                quiet = FALSE) {
  
  # create output directory and logfile
  dir.create(output_dir, showWarnings = FALSE)
  assign('logf', file(sprintf('%s/NMDer.log', output_dir), 'wt'), envir = .GlobalEnv)
  assign("quiet", quiet, envir = .GlobalEnv)
  options(warn=-1)
  
  # check for mandataory arguments
  if (any(missing(input), missing(reference), missing(fasta))) {
    stopLog('Missing mandatory arguments', logf)
  } 

  # import and/or load input file(s)
  c(inputGRanges, basicGRanges, genome) %<-% prepareInputs(input, reference, 
                                                           fasta, input_format, 
                                                           reference_format)

  # test for standard chromosome names, and gene_ids between input files
  inputGRanges = preTesting(inputGRanges, basicGRanges, genome, 
                            match_chrom, match_geneIDs, 
                            primary_gene_id, secondary_gene_id)
  
  # prepare dataframe and databases
  c(report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx) %<-% 
    prepareAnalysis(inputGRanges, basicGRanges, output_dir)
  
  # test each transcript for NMD features
  report_df = testNMDfeatures(report_df, inputExonsbyTx, basicExonsbyCDS,
                              basicExonsbyTx,genome, PTC_dist,other_features)
  
  # prepare outputs
  outputAnalysis(report_df, filterbycoverage, other_features, make_gtf, 
                 output_dir, input, reference, PTC_dist)
}

