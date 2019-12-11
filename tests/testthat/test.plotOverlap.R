

if (F){
  setwd('~/Toolbox/packages/pRoteomics/')  
  .libPaths(c("~/Toolbox/rlib/",.libPaths()))
  devtools::load_all()
  #library(testthat)
  library(ggplot2)
  library(gg)
}


test_that('test that overlap can be plotted',{
  
  # load data
  fpath = 'tests/testthat/data/14martin_FLT1_proteins.csv'
  dat = prepare(c('SMC','FLT'), fpath)
  colnames(dat)[1] <- 'gene'
  
  # for now sample from a genelist
  # find a bait db[db$gene %in% df$gene,]
  bait <- 'FLNA'
  db <- interactors(bait)
  db$significant <- FALSE
  plt = dat %>% mttest() %>% designate(FDR<1, logFC>0) %>% plotOverlap(bait, db, drawLabel = T, title = 'testing that FLNA can be drawn')
  plt
  expect_true(!is.null(plt))
  
})


test_that('test that Volcano plots can be plotted with multiple signifcant inweb interactors',{
  
  # load data
  fpath = 'tests/testthat/data/14martin_FLT1_proteins.csv'
  dat = prepare(c('SMC','FLT'), fpath)
  colnames(dat)[1] <- 'gene'
  
  # for now sample from a genelist
  # take out 10 random interactors that are inweb and the data (artifcially)
  # so that they can be plotted
  bait <- 'FLNA'
  db <- interactors(bait) # too many interactors.. let's subset them

  ## set threshold so that at least one interactor is significant
  plt = dat %>% mttest() %>% designate(FDR<1, logFC>1) %>% plotOverlap(bait, db, drawLabel = T)
  plt
  expect_true(!is.null(plt))
  
})



