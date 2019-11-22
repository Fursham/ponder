

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
    
    # Prepare list of reference transcript and GRanges
    thisline$Ref_transcript_ID = strsplit(thisline$Ref_transcript_ID, split = '_')
    queryGRanges = inputExonsbyTx %>% 
      dplyr::filter(group_name == thisline$Transcript_ID) %>%
      dplyr::arrange(ifelse(strand == '+', start, desc(start))) %>% 
      GenomicRanges::makeGRangesListFromDataFrame(keep.extra.columns = TRUE, split = 'group_name')
    
    basicTxGRanges = basicExonsbyTx %>% 
      dplyr::filter(group_name %in% thisline$Ref_transcript_ID[[1]]) %>%
      dplyr::arrange(ifelse(strand == '+', start, desc(start))) %>% 
      GenomicRanges::makeGRangesListFromDataFrame(keep.extra.columns = TRUE, split = 'group_name')

    # this function attempts to select the best reference for analysis 
    outBestRef = getBestRef(queryGRanges, basicTxGRanges)
    
    if(is.na(outBestRef$Ref_transcript_ID)){
      # return as query and reference do not match
      return(thisline)
    } else {
      # create new GRanges
      basicCDSGRanges = basicExonsbyCDS %>% 
        dplyr::filter(group_name %in% outBestRef$Ref_transcript_ID) %>%
        dplyr::arrange(ifelse(strand == '+', start, desc(start))) %>% 
        GenomicRanges::makeGRangesListFromDataFrame(keep.extra.columns = TRUE, split = 'group_name')
      basicTxGRanges = basicExonsbyTx %>% 
        dplyr::filter(group_name %in% outBestRef$Ref_transcript_ID) %>%
        dplyr::arrange(ifelse(strand == '+', start, desc(start))) %>% 
        GenomicRanges::makeGRangesListFromDataFrame(keep.extra.columns = TRUE, split = 'group_name')
      
      # update reference transcript and coverage
      thisline$Ref_transcript_ID = outBestRef$Ref_transcript_ID
      thisline$Coverage = outBestRef$Coverage
      
      # set query ORF if it is similar to reference
      if(outBestRef$Coverage == 1) {
        thisline$ORF_considered = basicCDSGRanges
        thisline$ORF_start = 'Annotated'
        thisline$ORF_found = TRUE
      }
     
    }
    
   
    
    # attempt to build Open Reading Frame for query
    ORFreport = getORF(basicCDSGRanges, queryGRanges,
                       genome, thisline$Gene_ID,
                       thisline$NMDer_ID)
    thisline = utils::modifyList(thisline, ORFreport)
    
    
    
    # if requested, test for NMD features and update line entry
    if (testforNMD == TRUE) {
      NMDreport = testNMD(thisline$ORF_considered, 
                          queryGRanges, 
                          PTC_dist, 
                          testNonClassicalNMD,
                          genome)
      thisline = utils::modifyList(thisline, NMDreport)
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
      
      thisline = utils::modifyList(thisline, altevents)
      
    }
    
    
    # update analyzed ORF coordinates into output
    #  this part have to be done this way because some terminal exons have length of 1
    #  and doing the conventional way will output 'xxx' rather than 'xxx-xxx'
    if (!is.na(thisline$ORF_considered[1])) {
    
      thisline$ORF_considered = thisline$ORF_considered %>% as.data.frame()
    }
    rm(list = c('queryGRanges','basicCDSGRanges', 'basicTxGRanges'))
    return(thisline)

  }
  
  # run the above function on report_df
  report_df = report_df %>% 
    dplyr::rowwise() %>% 
    do(data.frame(internalfunc(.), stringsAsFactors = FALSE)) 
  
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
  output = utils::modifyList(output, pre_report["annotatedStart"])
  
  # return if there is no shared exons between transcript and CDS
  if (is.na(pre_report$txrevise_out[1])) {
    return(output)
  } 
  
  # attempt to reconstruct CDS for transcripts with unannotated start
  if ((pre_report$annotatedStart == FALSE) |
      (pre_report$annotatedStart == TRUE & pre_report$firstexlen < 3)) {
    
    diffSegments = pre_report$txrevise_out
    queryStrand = as.character(strand(knownCDS))[1]
    pre_report$ORF = sort(append(
      diffSegments$shared_exons, c(
        reduce(diffSegments[[1]][diffSegments[[1]]$upstream != TRUE]),
        reduce(diffSegments[[2]][diffSegments[[2]]$upstream == TRUE]))),
      decreasing = queryStrand == '-')
    combinedList = list(refTx = pre_report$ORF, testTx = queryTx)
    pre_report$txrevise_out =  indentifyAddedRemovedRegions("refTx", "testTx", combinedList[c("refTx", "testTx")])
    
    pre_report$predictedStart = 'reconstructed'
    
    #pre_report = reconstructCDSstart(queryTx, knownCDS,
    #                                 refsequence,
    #                                 pre_report$txrevise_out,
    #                                 full.output = TRUE)
    output = utils::modifyList(output, list(predictedStart = pre_report$predictedStart))
    
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
  output = utils::modifyList(output, augmentedCDS)
  
  return(output)
}



getASevents <- function(transcript1, transcript2, testedNMD, orf, is_NMD) {
  
  # prepare list to be returned
  ASlist = data_frame(CE = NA, MX = NA, A5 = NA, A3 = NA, AF = NA, ATS = NA, AL = NA, APA = NA, IR = NA,
                      ce = NA, mx = NA, a5 = NA, a3 = NA, af = NA, ats = NA, al = NA, apa = NA, ir = NA) %>%
    dplyr::mutate_all(funs(as.character(.))) %>%
    as.list()
  if (testedNMD == TRUE) {
    ASlist = utils::modifyList(ASlist, list(NMDcausing = as.character(NA), NMDcausing.coord = as.character(NA)))
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
      ASlist = utils::modifyList(ASlist, list(NMDcausing = NMDexon, NMDcausing.coord = NMDexon.coord))
    }
  } 
  
  elementMetadata(combinedASoutput)$val = as.character(paste(ranges(combinedASoutput)))
  
  prepout = elementMetadata(combinedASoutput) %>% as.data.frame() %>%
    dplyr::group_by(AS_class) %>% 
    dplyr::summarise(vals = as.character(paste(val, collapse=";"))) %>% 
    as.data.frame()
  prepout2 = dplyr::select(prepout, vals) %>% unlist() %>% setNames(prepout[,1]) %>% as.list()
  ASlist = utils::modifyList(ASlist, prepout2)
  
  out = c(Shared_coverage = as.numeric(ASoutput$Shared_coverage), ASlist)
  return(out)
}




