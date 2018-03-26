context("Test txrevise")

test_that("identifyAddedRemoveRegions works", {
  diff = indentifyAddedRemovedRegions('ENSMUST00000029780','ENSMUST00000197833', ptbp2_testData$exons)
  
  expect_equal(length(diff$ENSMUST00000029780),3) #Has three elements
  expect_equal(width(diff$ENSMUST00000029780[diff$ENSMUST00000029780$contained == 1]), 34) #Width of the middle exon is 34 nt
})

context("Test NMDer")
test_that("testTxforNMD works", {
  NMDreport_ptbp2_noNMD = testTxforNMD(ptbp2_testData$noNMD, mmus_dna)
  NMDreport_ptbp2_NMD = testTxforNMD(ptbp2_testData$NMD, mmus_dna)
  NMDreport_bak1_noNMD = testTxforNMD(bak1_testData$noNMD, mmus_dna)
  NMDreport_bak1_NMD = testTxforNMD(bak1_testData$NMD, mmus_dna)
  
  
  # we need to store sequence info for Ptbp2 and Bak1 transcripts
  # i guess we can test multiple transcripts in this test function and have more expected output below
  
  expect_equal(mcols(NMDreport_ptbp2_noNMD[1])$lastEJ_dist,133)
  expect_equal(length(NMDreport_ptbp2_noNMD),1)
  
  expect_equal(mcols(NMDreport_ptbp2_NMD[3])$lastEJ_dist,-361)
  expect_equal(length(NMDreport_ptbp2_NMD),8)
  
  expect_equal(mcols(NMDreport_bak1_noNMD[1])$lastEJ_dist,105)
  expect_equal(length(NMDreport_bak1_noNMD),1)
  
  expect_equal(mcols(NMDreport_bak1_NMD[10])$lastEJ_dist,-89)
  expect_equal(length(NMDreport_bak1_NMD),22)
})