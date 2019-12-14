#' Match seqnames of input GRanges to reference GRanges
#'
#' @param from GRanges object with seqnames to change
#' @param to GRanges object from which seqnames is referenced
#'
#' @return Corrected input GRanges
#' @export
#' @author Fursham Hamid
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb mapSeqlevels
#'
matchSeqLevels <- function(from, to, ...){
  
  suppressWarnings(
    if (any(!seqlevelsStyle(from) %in% seqlevelsStyle(to))) {
      newStyle <- mapSeqlevels(seqlevels(from), (seqlevelsStyle(to)[1]))
      newStyle = newStyle[!is.na(newStyle)]
      from <- renameSeqlevels(from, newStyle)
      
      if (any(!seqlevels(from)%in%seqlevels(to))) {
        seqlevels(from, pruning.mode = 'tidy') <- as.vector(newStyle)
      }
    }
  )
  return(from)
}


