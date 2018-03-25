context("Test txrevise")

test_that("identifyAddedRemoveRegions works", {
  diff = indentifyAddedRemovedRegions('ENSMUST00000029780','ENSMUST00000197833', ptbp2_data$exons)
  
  expect_equal(length(diff$ENSMUST00000029780),3) #Has three elements
  expect_equal(width(diff$ENSMUST00000029780[diff$ENSMUST00000029780$contained == 1]), 34) #Width of the middle exon is 34 nt
})

context("Test NMDer")
test_that("testTxforNMD works", {
  NMDreport_ptbp2 = testTxforNMD(ptbp2_testData$cdss$ENSMUST00000029780, mmus_dna)
  NMDreport_bak1 = testTxforNMD(bak1_testData$cdss$ENSMUST00000078691, mmus_dna)
  
  # we need to store sequence info for Ptbp2 and Bak1 transcripts
  # i guess we can test multiple transcripts in this test function and have more expected output below
  
  expect_equal(mcols(NMDreport_ptbp2)$lastEJ_dist,133)
  expect_equal(mcols(NMDreport_bak1)$lastEJ_dist,105)
})