matchSeqLevels <- function(inGRanges, refGRanges){
  
  if (GenomeInfoDb::seqlevelsStyle(inGRanges) != GenomeInfoDb::seqlevelsStyle(refGRanges)) {
    newStyle <- GenomeInfoDb::mapSeqlevels(GenomeInfoDb::seqlevels(inGRanges), GenomeInfoDb::seqlevelsStyle(refGRanges))
    newStyle = newStyle[!is.na(newStyle)]
    inGRanges <- GenomeInfoDb::renameSeqlevels(inGRanges, newStyle)
    
    if (any(!GenomeInfoDb::seqlevels(inGRanges)%in%GenomeInfoDb::seqlevels(refGRanges))) {
      GenomeInfoDb::seqlevels(inGRanges, pruning.mode = 'tidy') <- as.vector(newStyle)
      warnLog('Non-standard chromosome IDs in query were removed')
    }
    
  }
  return(inGRanges)
}


