identifyNMDcausing <- function(ASgranges, cds){
  
  tx1index = c(1:length(ASgranges))
  tx2index = c((length(ASgranges)+1):(length(ASgranges)+length(cds)))
  strand = strand(ASgranges)[1] %>% as.character()
  

  disjoint = BiocGenerics::append(ASgranges, cds) %>%
    GenomicRanges::disjoin(with.revmap = T) %>%
    as.data.frame() %>%
    dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start)) %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(as = ifelse(any(revmap%in%tx1index), T,F)) %>%
    dplyr::mutate(cds = ifelse(any(revmap%in%tx2index), T,F)) %>%
    dplyr::mutate(stop = ifelse(any(revmap%in%tail(tx2index,1)), T,F)) %>%
    dplyr::mutate(olap = ifelse(lengths(revmap)==2, T,F)) %>%
    as.data.frame()
  disjoint = disjoint[(min(which(disjoint$cds == T))):nrow(disjoint),]
  
  disjoint = disjoint %>%
    dplyr::mutate(cdslead = dplyr::lead(cds, default = F)) %>%
    
  
  disjoint[(max(which(disjoint$type == 'cds'))+1):nrow(disjoint),]$type = 'down'
  
    
    dplyr::mutate(source = ifelse(revmap%in%tx1index, 1, 0)) %>%
    dplyr::mutate(type = ifelse(revmap%in%tx2index, 'olap','nolap')) %>%
    dplyr::mutate(source = ifelse(revmap%in%tx2index, 2, source))
  
  
    dplyr::mutate(type = ifelse(lengths(revmap)==2, 'olap','nolap')) 
  
  disjoint[1:min(which(disjoint$revmap %in% tx2index[1]))-1,]$type = 'up'
  disjoint[(max(which(disjoint$revmap %in% tail(tx2index,1)))+1):nrow(disjoint),]$type = 'down'
  
  %>%
    dplyr::mutate(source = ifelse(type != 'cons' & revmap%in%tx1index, 1, 0)) %>%
    dplyr::mutate(source = ifelse(type != 'cons' & revmap%in%tx2index, 2, source))
  
}