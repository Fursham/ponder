#' Predict sensitivity of mRNA transcripts to NMD
#'
#' @param tx 
#' GRanges object or GRangesList object containing exons
#' for each transcript. To search for NMD-inducing features, transcripts have
#' to be coding and thus contain a cds information of the same transcript name
#' @param cds 
#' GRanges object or GRangesList object containing coding regions (CDS)
#' for each transcript. GRangesList must have names that match names in tx,
#' else tx will not be analysed
#' @param NMDthreshold 
#' Minimum distance of STOP codon to last exon junction (EJ) which triggers NMD.
#' Default = 50bp
#' @param which 
#' List containing tx names to filter for analysis
#' @param return 
#' If tx and cds are GRangesList, returns only NMD-sensitive (default)
#' transcripts or all transcripts.
#'
#' @return
#' List with prediction of NMD sensitivity and statistics:
#' 
#' is_NMD: logical value in prediciting transcript sensitivity to NMD
#' 
#' dist_to_lastEJ: Integer value indicating distance of STOP codon to last EJ
#' Values are referenced from last EJ, thus a positive value indicates upstream
#' position of STOP codon while negative value indicates downstream position
#' 
#' num_of_down_EJs: Number of downstream EJs
#' 
#' dist_to_downEJs: Integer value indicating distance of STOP codon to each down 
#' EJs Values are referenced from last EJ, thus a positive value indicates 
#' upstream position of STOP codon while negative value indicates downstream 
#' position together with distances of STOP codon to last EJ
#' @export
#' @author Fursham Hamid
#'
#' @examples
#' 
#' ### To visualize transcripts
#' library(wiggleplotr)
#' plotTranscripts(query_exons, query_cds)
#' 
#' ### Examples with single GRanges objects
#' predictNMD(query_exons$transcript1, query_cds$transcript1) # NMD-insensitive
#' predictNMD(query_exons$transcript3, query_cds$transcript3) # NMD-sensitive
#' 
#' ### Examples with GRangesList object
#' predictNMD(query_exons, query_cds)
#' predictNMD(query_exons, query_cds), return = 'all')
#' predictNMD(query_exons, query_cds,which=c('transcript1', 'transcript3'), return = 'all')
#' 
#' 
#' 
predictNMD <- function(tx, cds, NMDthreshold = 50, 
                       which = NULL, return = c('NMD','all')){
  
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
  } else if(is(tx,'GRangesList') | is(tx,'list')){
    if(is(cds,'GRangesList') | is(cds,'list')){
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
      tx = tx[names(tx) %in% which]
    }
    # check for missing cds and return warnings/errors
    totest = totest[totest %in% names(cds)]
    if(length(totest)==0){
      stop('all tx have missing cds info. please ensure tx and cds names match')
    } 
    if(length(totest) < length(tx)){
      skiptest = length(tx) - length(totest)
      rlang::warn(sprintf('%s tx(s) have missing cds info and have been skipped',
                   skiptest))
    }
    
    # running brlapply and testNMD
    out = BiocParallel::bplapply(totest,  function(x){
      report = list(tx = x, is_NMD = F, dist_to_lastEJ = 0,
                    num_of_down_EJs = 0,dist_to_downEJs = 0)
      NMDreport = testNMD_(tx[[x]],cds[[x]],distance_stop_EJ=NMDthreshold)
      report = utils::modifyList(report, NMDreport)
      return(report)
    }, BPPARAM = BiocParallel::MulticoreParam()) %>%
      dplyr::bind_rows() %>%
      dplyr::filter(if(return[1]=='NMD') is_NMD ==T else T)
    return(out)
  }
}