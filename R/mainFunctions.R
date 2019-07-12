

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
    
    
    #############################################################################
    # deprecated
    
    # Add in coordinates of query transcripts
    #  this part have to be done this way because some exons are 1 nt in length
    #  and doing the conventional way will output 'xxx' rather than 'xxx-xxx'
    #string = paste(ranges(inputExonsbyTx[[thisline$Transcript_ID]]))
    #string = sapply(string, function(y) {
    #  newsstring = ifelse(!grepl('-', y), paste(c(y, '-', y), collapse = ''), y)
    #  return(newsstring)
    #})
    #thisline$Tx_coordinates = paste(string, collapse = ';')
    #############################################################################

    # this function attempts to select the best reference for analysis if multiple
    # reference transcripts exists
    if (length(thisline$Ref_TX_ID[[1]]) > 1) {
      # build GRanges object for query and GRangeslist object for reference
      queryGRanges = inputExonsbyTx %>% 
        filter(group_name == thisline$Transcript_ID) %>%
        makeGRangesListFromDataFrame(keep.extra.columns = TRUE, split = 'group_name')
      
      basicCDSGRanges = basicExonsbyCDS %>% 
        filter(group_name %in% thisline$Ref_TX_ID[[1]]) %>%
        makeGRangesListFromDataFrame(keep.extra.columns = TRUE, split = 'group_name')
      
      # calculate the coverage of query and reference transcripts as a percentage
      # of the mean width of query and reference
      overlapHits = mergeByOverlaps(queryGRanges, basicCDSGRanges)
      overlapHitsMeta = mapply(function(x,y){
        widthQuery = sum(width(x))
        widthRef = sum(width(y))
        aveWidth = (widthQuery + widthRef) / 2
        commonCoverage = sum(width(reduce(intersect(x, y))))
        Shared_coverage = commonCoverage/aveWidth
        return(Shared_coverage)
      }, overlapHits$queryGRanges, overlapHits$basicCDSGRanges) %>%
        as.data.frame() # output list as a dataframe
      
      # append reference tx IDs, sort dataframe and select best reference
      overlapHitsMeta$Ref_TX_ID = names(overlapHits$basicCDSGRanges)
      overlapHitsMeta = overlapHitsMeta %>% 
        filter(. < 1) %>%
        arrange(desc(.))
      
      thisline$Ref_TX_ID = overlapHitsMeta[1,2]
    } else {
      # convert list into character for comparisons with one reference transcript only
      thisline$Ref_TX_ID = as.character(thisline$Ref_TX_ID[[1]])
    }
    
    
    # prepare GRanges for analysis
    queryGRanges = inputExonsbyTx %>% 
      filter(group_name == thisline$Transcript_ID) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    basicCDSGRanges = basicExonsbyCDS %>% 
      filter(group_name == thisline$Ref_TX_ID) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    basicTxGRanges = basicExonsbyTx %>% 
      filter(group_name == thisline$Ref_TX_ID) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    # attempt to build Open Reading Frame for query
    ORFreport = getORF(basicCDSGRanges, queryGRanges,
                       genome, thisline$Gene_ID,
                       thisline$NMDer_ID)
    thisline = modifyList(thisline, ORFreport)
    
    
    # if requested, test for NMD features and update line entry
    if (testforNMD == TRUE) {
      NMDreport = testNMD(thisline$ORF_considered, 
                          queryGRanges, 
                          PTC_dist, 
                          testNonClassicalNMD,
                          genome)
      thisline = modifyList(thisline, NMDreport)
    }
    
    # if requested, classify alternative splicing events and update line entry
    if (testforAS == TRUE) {
      if (testforNMD == FALSE) {
        ORF = NA
        is_NMD = NA
      } else {
        ORF = thisline$ORF_considered
        is_NMD = thisline$is_NMD
      }
      
      altevents = getASevents(basicTxGRanges, queryGRanges, 
                              testforNMD, ORF, is_NMD)
      
      thisline = modifyList(thisline, altevents)
      
    }
    
    
    # update analyzed ORF coordinates into output
    #  this part have to be done this way because some terminal exons have length of 1
    #  and doing the conventional way will output 'xxx' rather than 'xxx-xxx'
    if (!is.na(thisline$ORF_considered[1])) {
      
      #### deprecated
      #ORFstringlist = paste(ranges(thisline$ORF_considered))
      #ORFstringlist = sapply(ORFstringlist, function(x){
      #  return(ifelse(!grepl('-', x), paste(c(x, '-', x), collapse = ''), x))
      #})
      #thisline$ORF_considered = paste(ORFstringlist, collapse = ';')
      ###################
      
      thisline$ORF_considered = thisline$ORF_considered %>% as.data.frame()
    }
    rm(list = c('queryGRanges','basicCDSGRanges', 'basicTxGRanges'))
    #thisline[11:length(thisline)] = as.character(thisline[11:length(thisline)])
    return(thisline)
    
    
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
getORF <- function(knownCDS, queryTx, refsequence, gene_id, transcript_id) {
  
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
  augmentedCDS = reconstructCDS(txrevise_out = pre_report$txrevise_out, 
                                fasta = refsequence, 
                                gene_id = gene_id, 
                                transcript_id = transcript_id)
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
    ASlist = modifyList(ASlist, list(NMDcausing = as.character(NA), NMDcausing.coord = as.character(NA)))
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
  strand = as.character(strand(transcript1)[1])
  combinedASoutput = append(ASoutput[[1]], ASoutput[[2]])
  combinedASoutput$AS_class = as.character(combinedASoutput$AS_class)
  combinedASoutput = sort(combinedASoutput, decreasing = strand == '-')
  combinedASoutput$size = width(combinedASoutput)
  
  NMDexon = NA
  NMDexon.coord = NA
  if (!is.na(is_NMD)) {
    if (is_NMD == TRUE) {
      
      # check if any of the alt segment overlaps with the last coding exon
      altseg_NMD = combinedASoutput[overlapsAny(combinedASoutput, orf[length(orf)])]
      if (length(altseg_NMD) == 1) {
        
        NMDexon = elementMetadata(altseg_NMD)$AS_class[1]
        NMDexon.coord = altseg_NMD[1] %>% range() %>% as.character() %>% substr(start = 1, stop = nchar(.)-2)
      } 
      # else, check if any of the alt segment overlaps with orf
      else {
        altseg_NMD = combinedASoutput[overlapsAny(combinedASoutput, range(orf))]
        
        # remove segments which are divisible by 3
        altseg_NMD = altseg_NMD[elementMetadata(altseg_NMD)$size %% 3 != 0]
        
        # if only 1 segment overlaps, that should be the NMD-causing exon
        if (length(altseg_NMD) == 1) {
          NMDexon = elementMetadata(altseg_NMD)$AS_class[1]
          NMDexon.coord = altseg_NMD[1] %>% range() %>% as.character() %>% substr(start = 1, stop = nchar(.)-2)
        } 
        # else, we need to test if each segment is in frame with orf
        else if (length(altseg_NMD) > 1) {
          NMDexon = paste(sort(elementMetadata(altseg_NMD)$AS_class), collapse = '|')
          NMDexon.coord = sort(altseg_NMD)[1] %>% range() %>% as.character() %>% substr(start = 1, stop = nchar(.)-2)
        
        }
        else if (length(altseg_NMD) == 0) {
          # if there is no internal out-of-frame exon causing the NMD, 
          # it could be due to spliced exons in 3'UTR that generate an exon-junction
          threeUTRseg = combinedASoutput[!overlapsAny(combinedASoutput, 
                                                      range(append(append(transcript1[1], transcript2[1]), 
                                                                   orf[length(orf)])))]
          
          if (length(threeUTRseg) > 0) {
            NMDexon = elementMetadata(threeUTRseg)$AS_class[1]
            NMDexon.coord = threeUTRseg[1] %>% range() %>% as.character() %>% substr(start = 1, stop = nchar(.)-2)
            NMDexon = paste(c('3UTR', NMDexon), collapse = '_')
          }
        }
      }
      ASlist = modifyList(ASlist, list(NMDcausing = NMDexon, NMDcausing.coord = NMDexon.coord))
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




