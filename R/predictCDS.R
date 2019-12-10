predictCDS <- function(query, refCDS, fasta,
                       query2ref, ids = c(1,2), 
                       coverage = NULL){
  
  # catch missing args
  mandargs <- ls()
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
  
  # create CDS list for tx with coverage of 1
  if(!is.null(coverage)){
    
  }
  
  
}
