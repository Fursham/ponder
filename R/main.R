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
run <- function(input,
                reference, 
                fasta,
                input_format = NULL,
                reference_format = NULL,
                match_chrom = FALSE,
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
  c(inputGRanges, basicGRanges, genome) %<-% prepareInputs(input, reference, fasta, input_format, reference_format)

  # test for standard chromosome names, and gene_ids between input files
  inputGRanges = preTesting(inputGRanges, basicGRanges, genome, correct_chrom, primary_gene_id, secondary_gene_id)
  
  # prepare dataframe and databases
  c(report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx) %<-% 
    prepareAnalysis(inputGRanges, basicGRanges, output_dir)
  
  # test each transcript for NMD features
  report_df = testNMDfeatures(report_df, 
                              inputExonsbyTx, basicExonsbyCDS,
                              basicExonsbyTx,
                              genome, PTC_dist,
                              other_features)
  

  # filter data based 
  if (filterbycoverage == TRUE) {
    report_df = filterdata(report_df)
  }
  
  # prepare splicing summary
  splicingsummarydf = summSplicing(report_df)
  write.table(splicingsummarydf, file = sprintf("%s/NMDer_splicing_summary.txt", output_dir), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  
  # prepare GTF file, if requested
  if (makegtf == TRUE) {
    generateGTF(report_df, output_dir)
  }
  
  # prepare output file
  infoLog('Saving analysis report...', logf, quiet)
  if (other_features == FALSE) {
    output_df = report_df[names(report_df) != c(uORF, threeUTR, uATG, uATG_frame)]
  } else {
    output_df = report_df
  }
  write(sprintf('# description; Input: %s; Reference: %s; PTC_to_EJ: %snt', 
                tail(unlist(strsplit(input, '/')), '1'), 
                tail(unlist(strsplit(reference, '/')), '1'), PTC_dist), 
        file = sprintf("%s/NMDer_report.txt", output_dir))
  write.table(output_df, file = sprintf("%s/NMDer_report.txt", output_dir), 
              sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
  unlink(logf)
}

