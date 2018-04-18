#' Run _ Program
#' 
#' @description
#' 
#'
#' @param file I
#' @param reference 
#' @param fasta 
#' @param primary_gene_id 
#' @param secondary_gene_id 
#' @param nonclassical 
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
run <- function(file,
                reference, 
                fasta,
                primary_gene_id = NULL,
                secondary_gene_id = NULL, 
                nonclassical = FALSE,
                makegtf = FALSE,
                PTC_dist = 50,
                output_dir = "NMDer",
                quiet = FALSE) {
  
  # check for mandataory arguments
  if (any(missing(file), missing(reference), missing(fasta))) {
    stop('Missing mandatory files and arguments')
  } 
  
  # create output directory and logfile
  dir.create(output_dir, showWarnings = FALSE)
  assign('logf', file(sprintf('%s/NMDer.log', output_dir), 'wt'), envir = .GlobalEnv)
  assign("quiet", quiet, envir = .GlobalEnv)
  options(warn=-1)

  # import and/or load input file(s)
  c(inputGRanges, basicGRanges, genome) %<-% prepareInputs(file, reference, fasta)

  # test for standard chromosome names, and gene_ids between input files
  inputGRanges = preTesting(inputGRanges, basicGRanges, genome, primary_gene_id, secondary_gene_id)
  
  # prepare output dataframe and databases
  c(report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx) %<-% 
    prepareAnalysis(inputGRanges, basicGRanges, output_dir)
  
  # test each transcript for NMD features
  report_df = testNMDfeatures(report_df, 
                              inputExonsbyTx, basicExonsbyCDS,
                              basicExonsbyTx,
                              genome, PTC_dist,
                              nonclassical)
  
  # prepare output file
  infoLog('Saving analysis report...', logf, quiet)
  #if (nonclassical == FALSE) {
  #  output_df = dplyr::select(report_df, Gene_ID, Ref_TX_ID, Original_Gene_ID, Gene_Name, 
  #                            ID_corrected, Transcript_ID, Chrom, Strand,
  #                            Tx_coordinates, annotatedStart, predictedStart, Alt_tx,
  #                            ORF_considered, is_NMD, dist_to_lastEJ)
  #} else {
  #  output_df = dplyr::select(report_df, Gene_ID, Ref_TX_ID, Original_Gene_ID, Gene_Name, 
  #                            ID_corrected, Transcript_ID, Chrom, Strand,
  #                            Tx_coordinates, annotatedStart, predictedStart, Alt_tx,
  #                            ORF_considered, is_NMD, dist_to_lastEJ, uORF, threeUTR)
  #}
  output_df = report_df
  
  write.table(output_df, file = sprintf("%s/NMDer_report.txt", output_dir), sep = "\t", row.names = FALSE)
  unlink(logf)
}

