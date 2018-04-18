context("Test txrevise")

test_that("identifyAddedRemoveRegions works", {
  diff = indentifyAddedRemovedRegions('ENSMUST00000029780','ENSMUST00000197833', ptbp2_testData$exons)
  
  expect_equal(length(diff$ENSMUST00000029780),3) #Has three elements
  expect_equal(width(diff$ENSMUST00000029780[diff$ENSMUST00000029780$contained == 1]), 34) #Width of the middle exon is 34 nt
})

context("Test testNMDvsThisCDS")
test_that("testNMDvsThisCDS", {
  ptbp2_out1 = testNMDvsThisCDS(ptbp2_testData$noNMD, ptbp2_testData$exons$ENSMUST00000197833, Ensembl_mm10)
  ptbp2_out2 = testNMDvsThisCDS(ptbp2_testData$noNMD, ptbp2_testData$exons$ENSMUST00000029780, Ensembl_mm10)
  ptbp2_out3 = testNMDvsThisCDS(ptbp2_testData$noNMD, ptbp2_testData$diffstart, Ensembl_mm10)
  
  expect_equal(ptbp2_out1$dist_to_lastEJ, -361)
  expect_equal(ptbp2_out2$Alt_tx, FALSE)
  expect_equal(ptbp2_out3$annotatedStart, FALSE)
})

context("Test testTXforStart")
test_that("testTXforStart", {
  out_ptbp2 = testTXforStart(ptbp2_testData$noNMD, ptbp2_testData$exons$ENSMUST00000197833, full.output=TRUE)
  out2_ptbp2 = testTXforStart(ptbp2_testData$noNMD, ptbp2_testData$diffstart, full.output=TRUE)
  
  expect_equal(length(out_ptbp2$txrevise_out$refTx), 1)
  expect_equal(out_ptbp2$annotatedStart, TRUE)
  expect_equal(length(out2_ptbp2$txrevise_out$refTx), 10)
  expect_equal(out2_ptbp2$annotatedStart, FALSE)
})

context("Test reconstructCDSstart")
test_that("reconstructCDSstart",  {
  newptbp2CDS = reconstructCDSstart(ptbp2_testData$noNMD, 
                                    ptbp2_testData$diffstart, 
                                    refsequence = Ensembl_mm10,
                                    full.output = TRUE)
  
  expect_equal(length(newptbp2CDS$txrevise_out$refTx), 9)
})



context("Test reconstructCDS")
test_that("reconstructCDS", {
  augmented_ptbp2 = reconstructCDS(ptbp2_testData$noNMD, ptbp2_testData$exons$ENSMUST00000197833)
  augmented_bak1 = reconstructCDS(bak1_testData$noNMD, bak1_testData$NMD)
  augmented_negative = reconstructCDS(ptbp2_testData$noNMD, ptbp2_testData$exons$ENSMUST00000029780)
  
  expect_equal(length(augmented_ptbp2$ORF), 13)
  expect_equal(length(augmented_bak1$ORF), 6)
  expect_equal(augmented_negative$Alt_tx, FALSE)  # test should return FALSE as transcript is a ref CDS
})


context("Test testClassicalNMD")
test_that("testClassicalNMD", {
  NMDreport_ptbp2_noNMD = testClassicalNMD(ptbp2_testData$noNMD, Ensembl_mm10)
  NMDreport_ptbp2_NMD = testClassicalNMD(ptbp2_testData$NMD, Ensembl_mm10)
  NMDreport_psd95_noNMD = testClassicalNMD(psd95_testData$noNMD, Ensembl_mm10)
  NMDreport_psd95_NMD = testClassicalNMD(psd95_testData$NMD, Ensembl_mm10)
  
  expect_equal(mcols(NMDreport_ptbp2_noNMD$stopcodons$lastEJ_dist),133)
  
  expect_equal(mcols(NMDreport_ptbp2_NMD$stopcodons[1])$lastEJ_dist,-361)
  expect_equal(length(NMDreport_ptbp2_NMD$stopcodons),8)
  
  expect_equal(length(NMDreport_psd95_noNMD$stopcodons),1)
  
  expect_equal(mcols(NMDreport_psd95_NMD$stopcodons[1])$lastEJ_dist,-80)
  expect_equal(length(NMDreport_psd95_NMD$stopcodons),2)
})

context("Test resizeTranscripts")
test_that("resizeTranscripts", {
  output = resizeTranscripts(ptbp2_testData$exons$ENSMUST00000197833, start = 100)
  
  expect_equal(length(output),11)
  expect_equal(end(output[1]),119761742)
})