#' Reconstruct CDS with alternative internal and downstream segments
#'
#' @description 
#' This function will add alternative segments from a query transcript into 
#' a reference CDS from the same gene and generate a new ORF
#' 
#' Note: This function will not insert/remove first exons. 
#' 
#' @param queryTranscript A GRanges object containing exon structure of a query transcript from
#' the same gene family
#' @param refCDS A GRanges object containing CCDS information of a gene family
#' @param fasta Fasta sequence of the genome
#' @param txrevise_out 
#' Optional. 
#' GRangesobject from 'indentifyAddedRemovedRegions()' 
#' or 'testTXforStart()' output. ?indentifyAddedRemovedRegions or ?testTXforStart
#' for more information.
#' Arguments refCDS and queryTranscript is not mandatory if this argument is provided
#' 
#' 
#' @return
#' A list containing:
#' (1) A GRanges Object of new ORF, or NA if no ORF is found
#' (2) TRUE/FALSE object on whether ORF is an alternative CDS transcript

#' @author Fursham Hamid
#'
#' @examples
#' 
#' library("BSgenome.Mmusculus.UCSC.mm10")
#' reconstructCDS(ptbp2Data$transcripts$ENSMUST00000197833, ptbp2Data$refCDS, fasta = BSgenome.Mmusculus.UCSC.mm10)
#' 
#' 
reconstructCDS <- function(queryTranscript, refCDS, fasta, txrevise_out = NULL, gene_id, transcript_id){
  
  # check if txrevise_out input is provided, and build one if not.
  if (is.null(txrevise_out)) {
    if (missing(refCDS) | missing(queryTranscript)) {
      stop('Please provide input GRanges objects')
    }
    combinedList = list(refTx = refCDS, testTx = queryTranscript)
    diffSegments = indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
  } else {
    diffSegments = txrevise_out
  }
  
  # prepare output
  output = list(ORF_considered = NA, ORF_found = FALSE)
  
  # return NA if there is no unique internal and downstream alternative segments
  # else, construct new CDS with insertion of segments from query transcript
  
  # Basically, this part tests if there is any unique internal and downstream segments
  if (length(diffSegments[[1]]$contained[diffSegments[[1]]$contained == TRUE]) == 0 &
      length(diffSegments[[1]]$downstream[diffSegments[[1]]$downstream == TRUE]) == 0 &
      length(diffSegments[[2]]$contained[diffSegments[[2]]$contained == TRUE]) == 0) {
    
    # if not, it will return the refCDS as the CDS and return alternative_tx as false
    Alternative_tx = FALSE
    augmentedCDS = sort(GenomicRanges::reduce(
      diffSegments$shared_exons),
      decreasing = as.character(strand(diffSegments$shared_exons))[1] == '-')
    
    ## POSSIBLE WARNING
  } else {
    # if there are, construct a new exon structure with removal/addition of the alternative segments 
    Alternative_tx = TRUE
    augmentedCDS = sort(append(
      diffSegments$shared_exons, 
      reduce(diffSegments[[2]][diffSegments[[2]]$upstream != TRUE])),
      decreasing = as.character(strand(diffSegments$shared_exons))[1] == '-')
    
    # this part will correct the open reading frame
    # obtain critical information on the new CDS
    queryStrand = as.character(strand(augmentedCDS))[1]
    thisqueryseq = unlist(Biostrings::getSeq(fasta, augmentedCDS))
    
    # prepare a dict of stop codons for pattern matching
    list_stopcodons = Biostrings::DNAStringSet(c("TAA", "TAG", "TGA"))
    pdict_stopcodons = Biostrings::PDict(list_stopcodons)
    
    # search for in-frame stop codons
    allmatches = Biostrings::matchPDict(pdict_stopcodons, thisqueryseq)
    combinedmatches = unlist(allmatches)
    inframe_stopcodons = sort(combinedmatches[end(combinedmatches) %% 3 == 0,])
    
    # append 3' end of transcript to the first stop codon if found
    if (length(inframe_stopcodons) > 0) {
      downUTRsize = length(thisqueryseq) - end(inframe_stopcodons[1])
      augmentedCDS = resizeTranscripts(augmentedCDS, end = downUTRsize)
      augmentedCDS = augmentedCDS %>% as.data.frame() %>%
        dplyr::mutate(type = 'CDS', gene_id = gene_id, transcript_id = transcript_id) %>%
        dplyr::mutate(phase = cumsum(width%%3)%%3)
      augmentedCDS$phase = c(0, head(augmentedCDS$phase, - 1))
      augmentedCDS = makeGRangesFromDataFrame(augmentedCDS, keep.extra.columns = TRUE)
      
      output$ORF_found = TRUE
    } else {
      # if a stop codon is not found, return NA
      augmentedCDS = NA
    }
  }
  output = modifyList(output, 
                      list(ORF_considered = augmentedCDS))
  return(output)
}

