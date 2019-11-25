context("Test prepNMDer workflow")

# perform workflow
library(BSgenome.Mmusculus.UCSC.mm10)
testQuerynew = paste0('../../', testQuery)
out = suppressWarnings(prepNMDer(testQuerynew, testRef, Mmusculus, match_chrom = T, match_geneIDs = T, primary_gene_id = 'gene_id', secondary_gene_id = 'ref_gene_id'))
test_that("Test the structure of output file", {
  
  expect_equal(isS4(out), TRUE) 
  expect_equal(length(slotNames(out)), 6)
  expect_equal(names(out@df), c('NMDer_ID', 'Gene_ID', 'Gene_name', 'Original_Gene_ID', 'Gene_match_level', 
                                'Transcript_ID', 'Transcript_coord', 
                                'Strand', 'Ref_transcript_ID'))
})

context("Test runNMDer workflow")

# perform workflow
out@df = out@df %>% dplyr::mutate(Coverage = 0,ORF_considered = as.character(NA), ORF_start = as.character('Not found'),ORF_found = FALSE)
out.main = suppressWarnings(runMain(out@df, out@inputTranscripts, out@basicCDS, out@basicTranscripts, Mmusculus, testforNMD = T))
test_that("Test the structure of output file", {
  
  expect_equal(names(out.main), names(runNMDer.out)) 
  expect_equal(out.main$Coverage, runNMDer.out$Coverage)
  expect_equal(out.main$is_NMD, runNMDer.out$is_NMD) 
})
