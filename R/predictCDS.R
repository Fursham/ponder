#' Title
#'
#' @param query 
#' @param refCDS 
#' @param fasta 
#' @param query2ref 
#' @param ids 
#' @param coverage 
#'
#' @return
#' @export
#'
#' @examples 
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' predictCDS(querytx, refCDS, Mmusculus, q2r_df)
#' 
predictCDS <- function(query, refCDS, fasta,
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
    querystyle = GenomeInfoDb::seqlevelsStyle(query)
    refstyle = GenomeInfoDb::seqlevelsStyle(refCDS)
    stop('query and refCDS has unmatched seqlevel styles. try matching using? function')
  }
  
  # sanity check if query and refCDS names are in q2f df
  if(all(!names(query) %in% query2ref[[ids[1]]])){
    unannotatedq = sum((!names(query) %in% query2ref[ids[1]]))
    rlang::warn(sprintf('%s query transcript ids were missing from query2ref df',
                        unnanotatedq))
  }
  if(all(!names(refCDS) %in% query2ref[[ids[2]]])){
    unannotatedr = sum((!names(refCDS) %in% query2ref[ids[2]]))
    rlang::warn(sprintf('%s reference CDS ids were missing from query2ref df',
                        unnanotatedr))
  }
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
    CDSreport = getORF(querytx[x], refCDS[y], fasta) %>%
      as.data.frame()
    return(CDSreport)
  },query2ref[[txname]], query2ref[[refname]],
  BPPARAM = BiocParallel::MulticoreParam()) %>%
    dplyr::bind_rows()
  outCDS = suppressWarnings(dplyr::bind_rows(outCDS, out) %>%
    dplyr::mutate(group_name = transcript_id) %>%
    GenomicRanges::makeGRangesListFromDataFrame(split.field = 'group_name',
                                                keep.extra.columns = T))
  
  return(outCDS)
}
