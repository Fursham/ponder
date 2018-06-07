

#' _ workflow: Test transcripts for NMD feature
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW.
#' This function will analyze the assembled transcripts for NMD features
#'
#' @param report_df 
#' dataframe of a list of assembled transcripts. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param basicTX_df 
#' dataframe of a list of reference CDS transcripts. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param inputExonsbyTx 
#' GRangesList of exon coordinates of assembled transcripts grouped by transcript name. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param basicExonsbyCDS 
#' GRangesList of exon coordinates of reference CDS grouped by transcript name. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param genome 
#' DNAstring containing genome sequence
#' @param PTC_dist 
#' Numerical value referring to minimium distance of premature stop codon to last exon junction
#' @param testNonClassicalNMD 
#' Boolean value (TRUE/FALSE) on whether to test for non-classical NMD features
#'
#' @return
#' updated form of report_df
#' @export
#' @import dplyr
#' @import GenomicRanges
#'
#' @examples
runMain <- function(report_df, inputExonsbyTx, basicExonsbyCDS, 
                            basicExonsbyTx, genome, 
                            testforNMD = TRUE, PTC_dist = 50, 
                            testNonClassicalNMD = FALSE, 
                            testforAS = FALSE) {
  
  # this is an internal function for testing NMD features on report_df rowwise
  internalfunc = function(x) {
    
    # convert each line into a list so that elements 
    # can be referenced as thisline$_
    thisline = as.list(x)
    
    # Add in coordinates of query transcripts
    #  this part have to be done this way because some exons are 1 nt in length
    #  and doing the conventional way will output 'xxx' rather than 'xxx-xxx'
    string = paste(ranges(inputExonsbyTx[[thisline$Transcript_ID]]))
    string = sapply(string, function(y) {
      newsstring = ifelse(!grepl('-', y), paste(c(y, '-', y), collapse = ''), y)
      return(newsstring)
    })
    thisline$Tx_coordinates = paste(string, collapse = ';')
    
    # return if query do not correspond to an annotated gene
    if (thisline$Match_level == 5) {
      return(thisline)
    } else {
      
      # return if reference transcript is not a CDS
      thisbasicCDS = basicExonsbyCDS[[thisline$Ref_TX_ID]]
      if (is.null(thisbasicCDS[1])) {
        return(thisline)
      } else {
        
        # run test and update output list
        ORFreport = getORF(thisbasicCDS, 
                           inputExonsbyTx[[thisline$Transcript_ID]], 
                           genome)
        thisline = modifyList(thisline, ORFreport)
        
        
        # test NMD only if requested
        if (testforNMD == TRUE) {
          NMDreport = testNMD(thisline$ORF_considered, 
                              inputExonsbyTx[[thisline$Transcript_ID]], 
                              PTC_dist, 
                              testNonClassicalNMD,
                              genome)
          thisline = modifyList(thisline, NMDreport)
        }
        
        # classify alternative splicing events if requested
        if (testforAS == TRUE) {
          if (testforNMD == FALSE) {
            ORF = NA
            is_NMD = NA
          } else {
            ORF = thisline$ORF_considered
            is_NMD = thisline$is_NMD
          }
          
          altevents = getASevents(basicExonsbyTx[[thisline$Ref_TX_ID]], 
                                  inputExonsbyTx[[thisline$Transcript_ID]], 
                                  testforNMD, ORF, is_NMD)
          
          thisline = modifyList(thisline, altevents)
          
        }
        
        
        # update analyzed ORF coordinates into output
        #  this part have to be done this way because some terminal exons have length of 1
        #  and doing the conventional way will output 'xxx' rather than 'xxx-xxx'
        if (!is.na(thisline$ORF_considered[1])) {
          ORFstringlist = paste(ranges(thisline$ORF_considered))
          ORFstringlist = sapply(ORFstringlist, function(x){
            return(ifelse(!grepl('-', x), paste(c(x, '-', x), collapse = ''), x))
          })
          thisline$ORF_considered = paste(ORFstringlist, collapse = ';')
        }
        
        #thisline[11:length(thisline)] = as.character(thisline[11:length(thisline)])
        return(thisline)
      }
    }
  }
  
  # run the above function on report_df
  report_df = report_df %>% 
    rowwise() %>% do(data.frame(internalfunc(.), stringsAsFactors = FALSE)) 
  
  return(report_df)
}




#' _ workflow: Test transcript for NMD feature against a reference CDS
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW.
#' This function will compare the query transcript with a reference CDS and tests
#' if insertion of alternative segments into the CDS generates an NMD substrate
#'
#' @param knownCDS 
#' @param queryTx 
#' @param refsequence 
#' @param PTC_dist 
#' @param nonClassicalNMD 
#'
#' @return
#' @export
#'
#' @examples
getORF <- function(knownCDS, queryTx, refsequence) {
  
  # prep output list
  output = list(ORF_considered = as.character(NA),
                Alt_tx  = as.logical(NA),
                annotatedStart = as.logical(NA),
                predictedStart = as.logical(NA))
  
  # precheck for annotated start codon on query transcript and update output
  pre_report = testTXforStart(queryTx, knownCDS, full.output=TRUE)
  output = modifyList(output, pre_report["annotatedStart"])
  
  # return if there is no shared exons between transcript and CDS
  if (is.na(pre_report$txrevise_out[1])) {
    return(output)
  } 
  
  # attempt to reconstruct CDS for transcripts with unannotated start
  if ((pre_report$annotatedStart == FALSE) |
      (pre_report$annotatedStart == TRUE & pre_report$firstexlen < 3)) {
    pre_report = reconstructCDSstart(queryTx, knownCDS,
                                     refsequence,
                                     pre_report$txrevise_out,
                                     full.output = TRUE)
    output = modifyList(output, list(predictedStart = pre_report$predictedStart))
    
    # return if CDS with new 5' do not contain a start codon
    if (is.na(pre_report$ORF[1])) {
      return(output)
    }
  } 
  
  # reconstruct CDS with insertion of alternative segments
  augmentedCDS = reconstructCDS(txrevise_out = pre_report$txrevise_out, fasta = refsequence)
  output = modifyList(output, augmentedCDS)
  
  return(output)
}



getASevents <- function(transcript1, transcript2, testedNMD, orf, is_NMD) {
  
  # prepare list to be returned
  ASlist = data_frame(CE = NA, MX = NA, A5 = NA, A3 = NA, AF = NA, ATS = NA, AL = NA, APA = NA, IR = NA,
                      ce = NA, mx = NA, a5 = NA, a3 = NA, af = NA, ats = NA, al = NA, apa = NA, ir = NA) %>%
    dplyr::mutate_all(funs(as.character(.))) %>%
    as.list()
  if (testedNMD == TRUE) {
    ASlist = modifyList(ASlist, list(NMDcausing = as.character(NA)))
  }
  
  # get AS classifications. transcript 1 is reference and transcript 2 is query in this case
  ASoutput = classifyAltSegments(transcript1, transcript2)
  
  # return if there is no alternative segments between transcripts
  if (is.null(ASoutput)){
    out = c(Shared_coverage = as.numeric(NA), ASlist)
    return(out)
  } else if (length(ASoutput[[1]]) == 0 & length(ASoutput[[2]]) == 0) {
    out = c(Shared_coverage = as.numeric(NA), ASlist)
    return(out)
  }
  
  # combine alternative segments and update ASlist with segment coordinates by matching class annotations with the named ASlist
  combinedASoutput = append(ASoutput[[1]], ASoutput[[2]])
  combinedASoutput$AS_class = as.character(combinedASoutput$AS_class)
  
  NMDexon = NA
  if (!is.na(is_NMD)) {
    if (is_NMD == TRUE) {
      
      # check if any of the alt segment overlaps with the last coding exon
      altseg_NMD = combinedASoutput[overlapsAny(combinedASoutput, orf[length(orf)])]
      if (length(altseg_NMD) == 1) {
        
        NMDexon = elementMetadata(altseg_NMD)$AS_class[1]
      } 
      # else, check if any of the alt segment overlaps with orf
      else {
        altseg_NMD = combinedASoutput[overlapsAny(combinedASoutput, range(orf))]
        
        # if only 1 segment overlaps, that should be the NMD-causing exon
        if (length(altseg_NMD) == 1) {
          NMDexon = elementMetadata(altseg_NMD)$AS_class[1]
        } 
        # else, we need to test if each segment is in frame with orf
        else if (length(altseg_NMD) > 1) {
          elementMetadata(altseg_NMD)$size = as.data.frame(width(altseg_NMD))
          
          # filter for in-frame segments
          altseg_NMD[elementMetadata(altseg_NMD)$size %% 3 == 0]
          
          # take the upstream most segment as the NMD-causing segment
          if (length(altseg_NMD) > 0) {
            NMDexon = elementMetadata(altseg_NMD)$AS_class[1]
          } 
        }
      }
      ASlist = modifyList(ASlist, list(NMDcausing = NMDexon))
    }
  } 
  
  elementMetadata(combinedASoutput)$val = as.character(paste(ranges(combinedASoutput)))
  
  prepout = elementMetadata(combinedASoutput) %>% as.data.frame() %>%
    dplyr::group_by(AS_class) %>% dplyr::summarise(vals = as.character(paste(val, collapse=";"))) %>% as.data.frame()
  prepout2 = dplyr::select(prepout, vals) %>% unlist() %>% setNames(prepout[,1]) %>% as.list()
  ASlist = modifyList(ASlist, prepout2)
  
  out = c(Shared_coverage = as.numeric(ASoutput$Shared_coverage), ASlist)
  return(out)
}




