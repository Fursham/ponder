context("Test predictNMD functions")

out1 = predictNMD(txl$ENSMUST00000029780,cdsl$ENSMUST00000029780)
out2 = predictNMD(txl,cdsl, return = 'all')

test_that("Test the structure of outputs", {
  
  expect_equal(names(out1), c('is_NMD', 'dist_to_lastEJ',
                              'num_of_down_EJs', 'dist_to_downEJs',
                              'threeUTRlength')) 
  expect_equal(names(out2), c('tx', 'is_NMD', 'dist_to_lastEJ',
                              'num_of_down_EJs', 'dist_to_downEJs',
                              'threeUTRlength')) 
  expect_equal(nrow(out2), 6) 

})

test_that("Test the value of outputs", {
  
  expect_equal(out2$dist_to_lastEJ, c(-133,361,-12,498,62,280)) 
  expect_equal(out2$num_of_down_EJs, c(0,3,0,1,1,3)) 
  expect_equal(out2$dist_to_downEJs, c("","66,283,361","",
                                       "498","62","56,232,280")) 
  expect_equal(out2$threeUTRlength, c(1586,502,343,736,788,381)) 

})
