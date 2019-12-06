getCoverage <- function(tx1, tx2){
  chrom = as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(tx1)))
  cov = GenomicRanges::coverage(c(tx1,tx2))
  index = which(names(cov) == chrom)
  cov = cov[[index]]
  cov_val = S4Vectors::runValue(cov)
  cov_len = S4Vectors::runLength(cov)
  return(sum(cov_len[cov_val==2]) / ((sum(cov_len[cov_val==1]) + (sum(cov_len[cov_val==2])*2))/2))
}