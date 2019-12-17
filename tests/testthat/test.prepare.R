# Install the released version from CRAN


if (F){
  ## manual testing
  setwd('~/Toolbox/packages/pRoteomics/')
  .libPaths(c("~/Toolbox/rlib/",.libPaths()))
  devtools::load_all()
  
}


## setup testing environment
fun = 'prepare'
data <- readRDS(file = 'tests/testthat/data/ms1.rds')


test_that('standard functionality',{
  
  # 
  id <- 'A1'
  res <- prepare(c('EC', 'EDNRA'), data)
  ref <- testpaths(id, fun)$ref
  expect_equal_to_reference(res, ref)
  #saveRDS(res, 'tests/testthat/reference/prepare/prepare.A1.rds')
  
})


## get targets official hgnc
uniprot <- acession.convert(acession.matrix(data$Accession)) 


test_that('proteins can be retained with filtering',{
  
  ## only one peptide
  id <- 'B1'
  target <- 'PEO1'
  target.uniprot <- uniprot[uniprot$symbol == target, ]$hgnc
  
  # with and without filtering
  wo <- c(target, target.uniprot) %in% prepare(c('EC', 'EDNRA'), data, firstcol = 'gene')$gene
  w <- c(target, target.uniprot) %in% prepare(c('EC', 'EDNRA'), data, firstcol = 'gene', filter.ignore = 'PEO1')$gene
  expect_false(any(wo)); expect_true(any(w))
  
  ## only one peptide
  id <- 'B2'
  target <- 'KAC4'
  target.uniprot <- uniprot[uniprot$symbol == target, ]$hgnc
  
  # with and without filtering
  wo <- c(target, target.uniprot) %in% prepare(c('EC', 'EDNRA'), data, firstcol = 'gene')$gene
  w <- c(target, target.uniprot) %in% prepare(c('EC', 'EDNRA'), data, firstcol = 'gene', filter.ignore = 'KAC4')$gene
  expect_false(any(wo)); expect_true(any(w))
  
})

test_that('Columns can be inputted manually',{
  
  
  # invalid cols generate an error
  id <- 'C1'
  cols = c('notvalidcolumn')
  expect_error(prepare(c('EC', 'EDNRA'), data, cols = cols))
  
  
  # Standard functionaliyu
  id <- 'C2'
  cols = c("Accession",
           'Intensity.EC.EDNRA1.iTRAQ4.114.', 'Intensity.EC.mock1.iTRAQ4.116.',
           'Intensity.EC.EDNRA2.iTRAQ4.115.', 'Intensity.EC.mock2.iTRAQ4.117.')
  res = prepare(c('EC', 'EDNRA'), data, cols = cols, verbose = T)
  ref = testpaths(id, fun)$ref
  expect_equal_to_reference(res, ref)
  
})






test_that('that this output matches something we have generated before',{
  
  ## check standard fucntionality
  id <- 'A1'
  bait = c('SMC','FLT')
  data = 'tests/testthat/data/14martin_FLT1_proteins.csv'
  dat = prepare(bait, data, impute = NULL,
                cols = c("Accession", "Intensity.SMC_FLT1.iTRAQ4.114.", "Intensity.SMC_mockF1.iTRAQ4.116.", 
                                     "Intensity.SMC_FLT2.iTRAQ4.115.", "Intensity.SMC_mockF2.iTRAQ4.117."))
  result <- dat[, c('rep1', 'rep2')]
  rownames(result) <- NULL
  
  reference <- read.csv('tests/testthat/reference/prepare/SMC_FLT1.NoImp.GenoppiInput.txt', sep = '\t')[,c(2,3)]
  reference <- reference[complete.cases(reference),]
  rownames(reference) <- NULL
  
  expect_equal(result,reference)
  
  ## check standard fucntionality while geussing columns in iTRAQ
  id <- 'A2'
  bait = c('SMC','FLT')
  data = 'tests/testthat/data/14martin_FLT1_proteins.csv'
  dat = prepare(bait, data)
  result <- dat[, c('rep1', 'rep2')]
  rownames(result) <- NULL
  
  reference <- read.csv('tests/testthat/reference/prepare/SMC_FLT1.NoImp.GenoppiInput.txt', sep = '\t')[,c(2,3)]
  reference <- reference[complete.cases(reference),]
  rownames(reference) <- NULL
  
  expect_equal(result,reference)
})

test_that('imputation is working: visual inspection',{
  
  ## check standard fucntionality
  id <- 'B1'
  bait = c('SMC','FLT')
  data = read.csv('tests/testthat/data/14martin_FLT1_proteins.csv')
  set.seed(1337)
  
  ## simulate some missing values
  s = sample(1:nrow(data), 250)
  data[s, as.vector(unlist(sapply(data, is.numeric)))] <- 0
  data$X.Unique <- 100
  
  ## 
  dat = prepare(bait, data, impute = list(stdwidth = 0.5, shift = -0.8), raw = T)

  x1 <- hist(dat[dat$imputed==1, ]$Intensity.SMC_FLT1.iTRAQ4.114., 100, xlim = c(-15, 15))
  x2 <- hist(dat[dat$imputed==0, ]$Intensity.SMC_FLT1.iTRAQ4.114., 100, xlim = c(-15, 15))
  
  valsWOImp = list(mu=mean(dat[dat$imputed==0, ]$rep1), sd=sd(dat[dat$imputed==0, ]$rep1))
  valsWImp = list(mu=mean(dat[dat$imputed==1, ]$rep1), sd=sd(dat[dat$imputed==1, ]$rep1))
  valsWOImp
  valsWImp
  
  # visual inspection
  plot(x2, col = 'grey', main = 'imputation') # Plot 1st histogram using a transparent color
  plot(x1, col = 'yellow', add = TRUE) # Add 2nd histogram using different color
  
})




