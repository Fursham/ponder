context("Test txrevise")

test_that("identifyAddedRemoveRegions works", {
  diff = indentifyAddedRemovedRegions('ENSMUST00000029780','ENSMUST00000197833', ptbp2_data$exons)
  
  expect_equal(length(diff$ENSMUST00000029780),3) #Has three elements
  expect_equal(width(diff$ENSMUST00000029780[diff$ENSMUST00000029780$contained == 1]), 34) #Width of the middle exon is 34 nt
})
