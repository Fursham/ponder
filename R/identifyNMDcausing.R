identifyNMDcausing <- function(ASgranges, cds){
  
  output = list(NMDcausing = NA, 
                NMDcausing.coord = NA)
  
  tx1index = c(1:length(ASgranges))
  tx2index = c((length(ASgranges)+1):(length(ASgranges)+length(cds)))
  lastcodingexonstart = start(cds[length(cds)])
  strand = strand(ASgranges)[1] %>% as.character()
  

  disjoint = BiocGenerics::append(ASgranges, cds) %>%
    GenomicRanges::disjoin(with.revmap = T) %>%
    as.data.frame() %>%
    dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start)) %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(as = ifelse(any(revmap%in%tx1index), T,F)) %>%
    dplyr::mutate(cds = ifelse(any(revmap%in%tx2index), T,F)) %>%
    dplyr::mutate(stop = ifelse(any(revmap%in%tail(tx2index,1)), 1,0)) %>%
    dplyr::mutate(outframe = ifelse(width%%3 > 0, 1,0)) %>%
    dplyr::mutate(dist = 0) %>%
    dplyr::mutate(append = '') %>%
    as.data.frame()
  disjoint = disjoint[(min(which(disjoint$cds == T))):nrow(disjoint),]
  disjoint[((max(which(disjoint$cds == T)))+1):nrow(disjoint),]$append = '3UTR-'
  disjoint$dist = rank(replace(disjoint$start - lastcodingexonstart, which(disjoint$start - lastcodingexonstart < 0), NA), na.last = T)
  
  
  disjoint = disjoint %>% dplyr::filter(as == T) %>%
    dplyr::arrange(dplyr::desc(stop), dplyr::desc(outframe), dist) %>%
    dplyr::left_join(ASgranges %>% as.data.frame(), by = c("seqnames", "start", "end", "width", "strand")) %>%
    dplyr::mutate(AS = paste0(append,AS)) %>%
    dplyr::select(seqnames:strand,AS)
  
  output$NMDcausing = disjoint[1,]$AS
  output$NMDcausing.coord = paste0(disjoint[1,2:3], collapse = '-')

  return(output)
}