testNMD <- function(tx, cds, NMDthreshold = 50, which = NULL){
  
  # catch missing args
  mandargs <- c('tx', 'cds')
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    stop(paste("missing values for", 
               paste(setdiff(mandargs, passed), collapse=", ")))
  }
  
  # check if tx and cds are GR or GRlist
  if(is(tx,'GRanges')){
    if(is(cds,'GRanges')){
      intype = 'gr'
    } else {
      txtype = is(tx)[1]
      cdstype = is(cds)[1]
      stop(sprintf('cds is type %s but tx is type %s',
                   cdstype, txtype))
    }
  } else if(is(tx,'list')){
    if(is(cds,'list')){
      intype = 'grl'
    } else {
      txtype = is(tx)[1]
      cdstype = is(cds)[1]
      stop(sprintf('cds is type %s but tx is type %s',
                   cdstype, txtype))
    }
  }
  
  #run testNMD_ for single GRanges object and output results
  if(intype=='gr'){
    return(testNMD_(tx,cds,distance_stop_EJ=NMDthreshold))
  }
  
  # for GRangesList, 
  if(intype=='grl'){
    totest = names(tx)  #prepare vector with names for testing
    if(!is.null(which)){
      totest = totest[totest %in% which] #subset list if which list is given
    }
    # check for missing cds and return warnings/errors
    skiptest = totest[!totest %in% names(cds)]
    if(length(skiptest)==length(totest)){
      stop('all tx have missing cds info. please ensure tx and cds names match')
    } else if(length(skiptest)>0){
      rlang::warn(sprintf('%s tx(s) have missing cds info and will be skipped',
                   length(skiptest)))
    }
    
    # running brlapply and testNMD
    out = BiocParallel::bplapply(totest,  function(x){
      report = list(tx = x, is_NMD = F, dist_to_lastEJ = 0)
      NMDreport = testNMD_(tx[[x]],cds[[x]],distance_stop_EJ=NMDthreshold)
      report = utils::modifyList(report, NMDreport)
      return(report)
    }, BPPARAM = BiocParallel::MulticoreParam()) %>%
      dplyr::bind_rows()
    return(out)
  }
}