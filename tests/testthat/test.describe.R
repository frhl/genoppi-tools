# Install the released version from CRAN

if (F){
  ## manual testing
  setwd('~/Toolbox/packages/pRoteomics/')
  .libPaths(c("~/Toolbox/rlib/",.libPaths()))
  devtools::load_all()
}

FUN = 'describe'
data <- readRDS(file = 'tests/testthat/data/ms1.rds')
test_that('iTRAQ columns are correctly identified',{
  
  set.seed(136)
  
  # standard functionality
  id <- 'A1'
  res <- describe(data)
  ref <- testpaths(id, FUN)$ref
  expect_equal_to_reference(res, ref)
  
  # scrambled columns are still captured correctly
  id <- 'A2'
  res <- data[,sample(ncol(data), 1:ncol(data))]
  ref <- testpaths(id, FUN)$ref
  expect_equal_to_reference(res, ref)
  
  # column names are all invalid
  id <- 'A3'
  res <- data
  colnames(res) <- 1:27
  res <- describe(res)
  expect_false(any(unlist(lapply(res, any))))
  
})





