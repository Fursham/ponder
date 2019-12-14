#' Construct query CDS using reference as guide
#'
#' @param query 
#' GRangesList object containing exons for each query transcript
#' @param refCDS 
#' GRangesList object containing CDS for each reference transcript
#' @param fasta 
#' BSgenome or Biostrings object containing genomic sequence
#' @param query2ref 
#' Dataframe with at least 2 columns: query transcript_id and its 
#' reference transcript_id. Query and ref transcript_ids have to match transcript
#' names in query and refCDS objects. IDs with missing GRanges object will
#' not be analysed
#' @param ids 
#' Numeric vector stating which columns of query2ref dataframe contain the 
#' query and reference transcript_ids respectively.
#' @param coverage 
#' Integer stating which column of query2ref dataframe contain percent coverage
#' between query and reference transcripts. Providing a column with coverage 
#' values will speed up CDS building process. Query transcripts that share 100%
#' coverage with reference CDS will be assigned the reference CDS and skip the
#' CDS searching process. See getCoverages function to calculate coverage values
#' (default:NULL)
#'
#' @return
#' GRangesList object containing CDS for each query transcript
#' @export
#'
#' @examples 
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' predictCDS(query, refCDS, Mmusculus, q2r_df)
#' 
buildCDS <- function(query, refCDS, fasta,
                       query2ref, ids = c(1,2), 
                       coverage = NULL){
  
  # catch missing args
  mandargs <- c('query','refCDS','fasta','query2ref')
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    stop(paste("missing values for", 
               paste(setdiff(mandargs, passed), collapse=", ")))
  }
  
  # catch unmatched seqlevels
  if(GenomeInfoDb::seqlevelsStyle(query) != GenomeInfoDb::seqlevelsStyle(refCDS)){
    stop('query and refCDS has unmatched seqlevel styles. try matching using matchSeqLevels function')
  }
  
  # sanity check if all tx in q2r have GRanges object
  if(!all(query2ref[[ids[1]]] %in% names(query))){
    missing = sum(!query2ref[[ids[1]]] %in% names(query))
    stop(sprintf('%s query transcripts have missing GRanges object',
                 missing))
  }
  if(!all(query2ref[[ids[2]]] %in% names(refCDS))){
    missing = sum(!query2ref[[ids[2]]] %in% names(refCDS))
    stop(sprintf('%s reference CDSs have missing GRanges object',
                 missing))
  }
  
  ### Not sure if the check below is necessary
  # # sanity check if query and ref names are in q2r df
  # if(all(!names(query) %in% query2ref[[ids[1]]])){
  #   unannotatedq = sum((!names(query) %in% query2ref[ids[1]]))
  #   rlang::warn(sprintf('%s query transcript ids were missing from query2ref df',
  #                       unannotatedq))
  # }
  # if(all(!names(refCDS) %in% query2ref[[ids[2]]])){
  #   unannotatedr = sum((!names(refCDS) %in% query2ref[ids[2]]))
  #   rlang::warn(sprintf('%s reference CDS ids were missing from query2ref df',
  #                       unannotatedr))
  # }
  
  # extract colnames and prepare outputCDS
  txname = names(query2ref)[ids[1]]
  refname = names(query2ref)[ids[2]]
  outCDS = data.frame()
  
  # create CDS list for tx with coverage of 1
  if(!is.null(coverage)){
    covname = names(query2ref)[coverage]  #extract cov colname
    #subset q2r for full coverages
    fullcovs = query2ref %>% 
      dplyr::filter(!!as.symbol(covname) == 1)
    query2ref = query2ref %>% 
      dplyr::filter(!!as.symbol(covname) != 1)
    # prepare outputCDS for full coverages
    outCDS = refCDS %>%
      as.data.frame() %>%
      dplyr::filter(group_name %in% fullcovs[[ids[2]]]) %>%
      dplyr::select(group_name:strand) %>%
      dplyr::mutate(type = 'CDS') %>%
      dplyr::left_join(fullcovs[ids], 
                       by = c('group_name'=refname)) %>%
      dplyr::mutate(phase = cumsum(width%%3)%%3) %>%
      dplyr::mutate(phase = dplyr::lag(phase, default = 0)) %>%
      dplyr::select(-group_name)
  }
  
  # create CDS list for all remaining tx
  out = BiocParallel::bpmapply(function(x,y){
    CDSreport = getORF(query[x], refCDS[y], fasta) %>%
      as.data.frame()
    return(CDSreport)
  },query2ref[[txname]], query2ref[[refname]],
  BPPARAM = BiocParallel::MulticoreParam()) %>%
    dplyr::bind_rows()
  outCDS = suppressWarnings(dplyr::bind_rows(outCDS, out) %>%
    dplyr::mutate(group_name = transcript_id) %>%
    GenomicRanges::makeGRangesListFromDataFrame(split.field = 'group_name',
                                                keep.extra.columns = T))
  
  # warn users if program fails to find CDS for some transcripts
  if(length(outCDS) < length(query)){
    missingCDS = length(query) - length(outCDS)
    rlang::warn(sprintf('Unable to find CDS for %s transcripts',missingCDS))
  }
  
  return(outCDS)
}
