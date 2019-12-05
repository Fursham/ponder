getCoverage <- function(tx1, tx2){
  chrom = as.character(runValue(seqnames(tx1)))
  cov = coverage(c(tx1,tx2))
  index = which(names(cov) == chrom)
  cov = cov[[index]]
  cov_val = runValue(cov)
  cov_len = runLength(cov)
  return(sum(cov_len[cov_val==2]) / ((sum(cov_len[cov_val==1]) + (sum(cov_len[cov_val==2])*2))/2))
}