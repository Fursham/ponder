getCoverages <- function(query, ref, query2ref, 
                         ids = c(1,2), return = c('best','all')){
  
  #Plans: ensure query2ref ids are in query and refCDS
  
  # catch missing args
  mandargs <- c('query','ref','query2ref')
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    stop(paste("missing values for", 
               paste(setdiff(mandargs, passed), collapse=", ")))
  }
  
  # catch unmatched seqlevels
  if(GenomeInfoDb::seqlevelsStyle(query) != GenomeInfoDb::seqlevelsStyle(ref)){
    querystyle = GenomeInfoDb::seqlevelsStyle(query)
    refstyle = GenomeInfoDb::seqlevelsStyle(ref)
    stop('query and ref has unmatched seqlevel styles. try matching using? function')
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
  
  # sanity check if query and ref names are in q2f df
  if(all(!names(query) %in% query2ref[[ids[1]]])){
    unannotatedq = sum((!names(query) %in% query2ref[ids[1]]))
    rlang::warn(sprintf('%s query transcript ids were missing from query2ref df',
                        unnanotatedq))
  }
  if(all(!names(ref) %in% query2ref[[ids[2]]])){
    unannotatedr = sum((!names(ref) %in% query2ref[ids[2]]))
    rlang::warn(sprintf('%s reference CDS ids were missing from query2ref df',
                        unnanotatedr))
  }
  
  # extract colnames and prepare outputCDS
  txname = names(query2ref)[ids[1]]
  refname = names(query2ref)[ids[2]]
  
  # get Coverage values for all comparisons
  out = BiocParallel::bpmapply(function(x,y){
    covrep = getCoverage(query[[x]], ref[[y]])
    return(covrep)
  },query2ref[[txname]], query2ref[[refname]],
  BPPARAM = BiocParallel::MulticoreParam())
  query2ref$coverage = out
  
  if(return[1] == 'best'){
    query2ref = query2ref %>%
      dplyr::arrange(!!as.symbol(txname), dplyr::desc(coverage)) %>%
      dplyr::distinct(!!as.symbol(txname), .keep_all = T)
  }
  return(query2ref)
}