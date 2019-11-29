
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
  ASoutput = classifyAS(transcript1, transcript2)
  strand = as.character(strand(transcript1)[1])
  

  # combine alternative segments and update ASlist with segment coordinates by matching class annotations with the named ASlist
  
  ASoutput = ASoutput %>% as.data.frame()

  
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

