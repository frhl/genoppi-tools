
if (F){
  ## manual testing
  setwd('~/Toolbox/packages/pRoteomics/')  
  .libPaths(c("~/Toolbox/rlib/",.libPaths()))
  devtools::load_all()
  library(testthat)
  
}



test_that('test that Volcano plots can be plotted',{
  
  # load data
  fpath = 'tests/testthat/data/14martin_FLT1_proteins.csv'
  dat = prepare(c('SMC','FLT'), fpath)
  colnames(dat)[1] <- 'gene'
  plt = dat %>% mttest() %>% designate(FDR<0.1) %>% plotVolcano('SMC')
  expect_true(!is.null(plt))
  
})


