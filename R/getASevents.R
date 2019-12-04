
getASevents <- function(transcript1, transcript2, testedNMD, orf, is_NMD) {
  
  # prepare list to be returned
  ASlist = list(CE = 'NA', MX = 'NA', A5 = 'NA', A3 = 'NA', AF = 'NA', ATS = 'NA', AL = 'NA', APA = 'NA', IR = 'NA',
                      ce = 'NA', mx = 'NA', a5 = 'NA', a3 = 'NA', af = 'NA', ats = 'NA', al = 'NA', apa = 'NA', ir = 'NA') 
  if (testedNMD == TRUE) {
    ASlist = utils::modifyList(ASlist, list(NMDcausing = as.character('NA'), NMDcausing.coord = as.character('NA')))
  }
  
  ASreport = classifyAS(transcript1, transcript2)
  ASreport_out = ASreport %>% as.data.frame() %>% 
    dplyr::group_by(AS) %>% 
    dplyr::mutate(coord = paste0(start,'-',end)) %>% 
    dplyr::summarise(combcoord = as.character(paste(coord, collapse = ';'))) %>% 
    dplyr::select(AS, combcoord)
  ASlist = utils::modifyList(ASlist, split(ASreport_out$combcoord, ASreport_out$AS))
  
  
  if(testedNMD == T & is_NMD == T){
    NMDexonreport = identifyNMDcausing(ASreport, orf)
    ASlist = utils::modifyList(ASlist, NMDexonreport)
  }
  
  return(ASlist)
}

