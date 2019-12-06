#' Resize GRanges transcript object
#' 
#' @description 
#' This function will append/extend the 5' and 3' ends of a GRanges object which describes the 
#' exon ranges of a transcript or a CDS
#'
#' @usage resizeTranscripts(x, start = 0, end = 0)
#'
#' @param x GRanges object containing exon coordinates of a transcript or CDS
#' @param start Length of 5' end to truncate (positive val) or extend (negative val)
#' @param end Length of 3' end to truncate (positive val) or extend (negative val)
#'
#' @return a new GRanges transcript object
#'
#' @examples
#' 
#' resizeTranscripts(ptbp2Data$transcripts$ENSMUST00000197833, start = 100)
#' 
old_resizeTranscripts <- function(x, start = 0, end = 0) {
  
  if (sum(width(x)) < (start + end)) {
    stop('Appending length is larger than size of transcript')
  }
  
  # retrieve important information
  strand = as.character(strand(x))[1]
  x = sort(x, decreasing = strand == '-')
  firstlastexons = c(1, length(x))  # initiate index of first and last exons
  startend = c(start,end) # store value for appending
  
  totest = sapply(startend, function(x) {
    ifelse(x > 0, TRUE, FALSE)
  })
  # append the start followed by the end
  for (i in (1:2)[totest]) {
    exonsize = width(x) # get the size of each exons as a list
    loop = TRUE
    while (loop == TRUE) {
      # change the size of the first or last exon
      exonsize[firstlastexons[[i]]] = exonsize[firstlastexons[[i]]] - startend[[i]]
      
      # if all of the exon sizes are more than 0, it means that the subtraction occurs within an exon
      #   if there are some exon with negative sizes, it means that the subtraction bleed into the next exon
      if (all(exonsize >= 0)) {
        loop = FALSE
      } else {
        
        # set the remaining start/end value to be subtracted
        newval = ifelse(exonsize[firstlastexons[[i]]] < 0, abs(exonsize[firstlastexons[[i]]]), 0)
        startend[[i]] = newval
        
        # convert negative widths to 0
        exonsize = ifelse(exonsize > 0, exonsize, 0)
        
        # and change index of exon to be subtracted
        if (i == 1) {
          firstlastexons[1] = ifelse(exonsize[firstlastexons[1]] > 0, firstlastexons[1], firstlastexons[1] + 1)
        } else if (i == 2) {
          firstlastexons[2] = ifelse(exonsize[firstlastexons[2]] > 0, firstlastexons[2], firstlastexons[2] - 1)
        }
      }
      
    }
    # resize the transcripts based on its new width
    if (i == 1) {
      x[1:length(x)] = resize(x[1:length(x)], 
                              exonsize[1:length(x)], 
                              fix = ifelse(strand == '-', "start", "end"), 
                              ignore.strand = TRUE)
    } else if (i == 2) {
      x[1:length(x)] = resize(x[1:length(x)], 
                              exonsize[1:length(x)], 
                              fix = ifelse(strand == '-', "end", "start"), 
                              ignore.strand = TRUE)
    }
  }
  # return exons with width of more than 0
  x = x[width(x) > 0]
  return(x)
}

