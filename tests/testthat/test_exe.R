context("Test prepNMDer workflow")

# perform workflow
out = prepNMDer(testQuery, testRef, Mmusculus, match_chrom = T, match_geneIDs = T, primary_gene_id = 'gene_id', secondary_gene_id = 'ref_gene_id')
test_that("Test the structure of output file", {
  
  expect_equal(isS4(out), TRUE) 
  expect_equal(length(slotNames(out)), 6)
  expect_equal(names(out@df), c('NMDer_ID', 'Gene_ID', 'Original_Gene_ID', 
                                'Match_level', 'Gene_name', 'Transcript_ID', 
                                'Ref_transcript_ID', 'Transcript_coord', 'Strand'))
})
