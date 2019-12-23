context("Test searchuORFs functions")

library(BSgenome.Mmusculus.UCSC.mm10)
out = searchuORFs(query_exons, query_cds, Mmusculus)

test_that("Test the structure of outputs", {
  
  expect_equal(names(out), c('exons', 'cds'))
  expect_equal(names(out$exons), c("uORF_1_transcript2", "uORF_2_transcript2"))
  expect_equal(sum(lengths(out$exons)), 28)
  expect_equal(sum(lengths(out$cds)), 2)
})

test_that("Test the value of outputs", {
  
  expect_equal(start(out$exons)[[1]][1:3], c(119783202, 119781504, 119761703)) 
  expect_equal(end(out$exons)[[1]][1:3], c(119783804, 119781534, 119761778)) 
  
  expect_equal(start(out$cds)[[1]], 119783744)
  expect_equal(end(out$cds)[[1]], 119783803)
  
  expect_equal(out$cds[[2]]$frame, 2)
  expect_equal(out$cds[[1]]$phase, 0)

})
