context("Test txrevise")

test_that("identifyAddedRemoveRegions works", {
  diff = indentifyAddedRemovedRegions('ENSMUST00000029780','ENSMUST00000197833', ptbp2Data$transcripts)
  
  expect_equal(length(diff$ENSMUST00000029780),3) #Has three elements
  expect_equal(width(diff$ENSMUST00000029780[diff$ENSMUST00000029780$contained == 1]), 34) #Width of the middle exon is 34 nt
})


context("Test testTXforStart")
test_that("testTXforStart", {
  out_ptbp2 = testTXforStart(ptbp2Data$transcripts$ENSMUST00000197833, ptbp2Data$refCDS, TRUE)
  out2_ptbp2 = testTXforStart(ptbp2Data$afCDS, ptbp2Data$refCDS, TRUE)
  
  #' testTXforStart(ptbp2Data$transcripts$ENSMUST00000197833, ptbp2Data$refCDS)
  #' testTXforStart(ptbp2Data$afCDS, ptbp2Data$refCDS)
  
  expect_equal(length(out_ptbp2$txrevise_out$refTx), 1)
  expect_equal(out_ptbp2$annotatedStart, TRUE)
  expect_equal(length(out2_ptbp2$txrevise_out$refTx), 10)
  expect_equal(out2_ptbp2$annotatedStart, FALSE)
})

context("Test reconstructCDSstart")
test_that("reconstructCDSstart",  {
  newptbp2CDS = reconstructCDSstart(ptbp2Data$afCDS, ptbp2Data$refCDS, fasta = BSgenome.Mmusculus.UCSC.mm10, full.output = TRUE)
  
  expect_equal(length(newptbp2CDS$txrevise_out$refTx), 9)
  expect_equal(newptbp2CDS$predictedStart, TRUE)
})



context("Test reconstructCDS")
test_that("reconstructCDS", {
  augmented_ptbp2 = reconstructCDS(ptbp2Data$transcripts$ENSMUST00000197833, ptbp2Data$refCDS, fasta = BSgenome.Mmusculus.UCSC.mm10)
  augmented_bak1 = reconstructCDS(bak1Data$transcripts$ENSMUST00000025034, bak1Data$refCDS, fasta = BSgenome.Mmusculus.UCSC.mm10)
  augmented_negative = reconstructCDS(ptbp2Data$transcripts$ENSMUST00000029780, ptbp2Data$refCDS, fasta = BSgenome.Mmusculus.UCSC.mm10)
  
  expect_equal(length(augmented_ptbp2$ORF_considered), 10)
  expect_equal(augmented_ptbp2$Alt_tx, TRUE)
  expect_equal(length(augmented_bak1$ORF), 5)
  expect_equal(augmented_bak1$Alt_tx, TRUE)
  expect_equal(augmented_negative$Alt_tx, FALSE)  # test should return FALSE as transcript is a ref CDS
})


context("Test testNMD")
test_that("testNMD", {
  NMDreport_ptbp2_noNMD = testNMD(ptbp2Data$refCDS, ptbp2Data$transcripts$ENSMUST00000029780)
  NMDreport_ptbp2_noNMDfull = testNMD(ptbp2Data$refCDS, ptbp2Data$transcripts$ENSMUST00000029780, other_features = TRUE, fasta = BSgenome.Mmusculus.UCSC.mm10)
  NMDreport_ptbp2_NMD = testNMD(ptbp2Data$skipE10CDS, ptbp2Data$transcripts$ENSMUST00000197833)
  
  NMDreport_psd95_noNMD = testNMD(psd95Data$refCDS, psd95Data$transcripts$ENSMUST00000108589)
  NMDreport_psd95_NMD = testNMD(psd95Data$poisonCDS, psd95Data$transcripts$ENSMUST00000123687)
  
  expect_equal(NMDreport_ptbp2_noNMD$dist_to_lastEJ,133)
  expect_equal(NMDreport_ptbp2_noNMDfull$threeUTR,1586)
  expect_equal(NMDreport_ptbp2_noNMDfull$uORF,FALSE)
  expect_equal(NMDreport_ptbp2_NMD$dist_to_lastEJ,-361)

  expect_equal(NMDreport_psd95_noNMD$is_NMD, FALSE)
  expect_equal(NMDreport_psd95_NMD$dist_to_lastEJ,-80)
})

context("Test resizeTranscripts")
test_that("resizeTranscripts", {
  output = resizeTranscripts(ptbp2Data$transcripts$ENSMUST00000197833, start = 100)
  
  expect_equal(length(output),11)
  expect_equal(end(output[1]),119761742)
})

context("Test classifyAltSegments")
test_that("classifyAltSegments", {
  output = classifyAltSegments(ptbp2Data$transcripts$ENSMUST00000029780, ptbp2Data$transcripts$ENSMUST00000197833)

  expect_equal(as.character(elementMetadata(output$transcript1)[2,]),'ce')

})