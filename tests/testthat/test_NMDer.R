context("Test txrevise")

test_that("identifyAddedRemoveRegions works", {
  diff = indentifyAddedRemovedRegions('ENSMUST00000029780','ENSMUST00000197833', ptbp2_testData$exons)
  
  expect_equal(length(diff$ENSMUST00000029780),3) #Has three elements
  expect_equal(width(diff$ENSMUST00000029780[diff$ENSMUST00000029780$contained == 1]), 34) #Width of the middle exon is 34 nt
})

context("Test NMDer")
test_that("augmentCDS", {
  augmented_ptbp2 = augmentCDS(ptbp2_testData$noNMD, exons$ENSMUST00000197833)
  augmented_bak1 = augmentCDS(bak1_testData$noNMD, exons$ENSMUST00000025034)
  augmented_negative = augmentCDS(ptbp2_testData$noNMD, exons$ENSMUST00000029780)
  
  expect_equal(length(augmented_ptbp2), 13)
  expect_equal(length(augmented_bak1), 6)
  expect_equal(augmented_negative, NULL)  # test should return NULL as transcript is a ref CDS
})
test_that("testTxforNMD", {
  NMDreport_ptbp2_noNMD = testTxforNMD(ptbp2_testData$noNMD, mmus_dna)
  NMDreport_ptbp2_NMD = testTxforNMD(ptbp2_testData$NMD, mmus_dna)
  NMDreport_psd95_noNMD = testTxforNMD(psd95_testData$noNMD, mmus_dna)
  NMDreport_psd95_NMD = testTxforNMD(psd95_testData$NMD, mmus_dna)
  
  
  # we need to store sequence info for Ptbp2 and Bak1 transcripts
  # i guess we can test multiple transcripts in this test function and have more expected output below
  
  expect_equal(mcols(NMDreport_ptbp2_noNMD[1])$lastEJ_dist,133)
  expect_equal(length(NMDreport_ptbp2_noNMD),1)
  
  expect_equal(mcols(NMDreport_ptbp2_NMD[3])$lastEJ_dist,-361)
  expect_equal(length(NMDreport_ptbp2_NMD),8)
  
  expect_equal(mcols(NMDreport_psd95_noNMD[1])$lastEJ_dist,107)
  expect_equal(length(NMDreport_psd95_noNMD),1)
  
  expect_equal(mcols(NMDreport_psd95_NMD[1])$lastEJ_dist,-80)
  expect_equal(length(NMDreport_psd95_NMD),2)
})