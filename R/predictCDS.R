predictCDS <- function(query, ref, fasta,
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
  if(GenomeInfoDb::seqlevelsStyle(query) != GenomeInfoDb::seqlevelsStyle(ref)){
    querystyle = GenomeInfoDb::seqlevelsStyle(query)
    refstyle = GenomeInfoDb::seqlevelsStyle(ref)
    stop('query and ref has unmatched seqlevel styles. try matching using? function')
  }
  
  # sanity check if query and ref names are in q2f df
  if(all(!names(query) %in% query2ref[ids[1]])){
    unannotatedq = sum((!names(query) %in% query2ref[ids[1]]))
    rlang::warn(sprintf('%s query transcripts were missing from query2ref df',
                        unnanotated))
  }
  if(all(!names(ref) %in% query2ref[ids[2]])){
    unannotatedq = sum((!names(ref) %in% query2ref[ids[2]]))
    rlang::warn(sprintf('%s ref CDS were missing from query2ref df',
                        unnanotated))
  }
  
  # create CDS list for tx with coverage of 1
  if(!is.null(coverage)){
    
  }
  
  
}
