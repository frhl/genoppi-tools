# Install the released version from CRAN

.libPaths(c("~/Toolbox/rlib/",.libPaths()))
devtools::load_all()



test_that('that this output matches something we have generated before',{
  
  ## check standard fucntionality
  id <- 'A1'
  bait = c('SMC','FLT')
  data = 'tests/testthat/data/14martin_FLT1_proteins.csv'
  dat = prepare(bait, data, cols = c("Accession", "Intensity.SMC_FLT1.iTRAQ4.114.", "Intensity.SMC_mockF1.iTRAQ4.116.", 
                                     "Intensity.SMC_FLT2.iTRAQ4.115.", "Intensity.SMC_mockF2.iTRAQ4.117."))
  result <- dat[, c('rep1', 'rep2')]
  reference <- read.csv('tests/testthat/reference/prepare/SMC_FLT1.NoImp.GenoppiInput.txt', sep = '\t')[,c(2,3)]
  expect_equal(result,reference)
})


# need to do it bait-wise instead of file-wise
#bait = unlist(strsplit(c('SMC_EDNRA','SMC_FLT','EC_FLT1')[2],'\\_')) #,'EC_JCAD')
#infile = list.files('~/Desktop/MICOM/WhiteheadData/', full.names = T)[2]
#dat = prepare(bait, infile, cols = c("Accession", "Intensity.SMC_FLT1.iTRAQ4.114.", "Intensity.SMC_mockF1.iTRAQ4.116.", 
#                                            "Intensity.SMC_FLT2.iTRAQ4.115.", "Intensity.SMC_mockF2.iTRAQ4.117."))

# 0.45441764  0.10679001  0.89405064 -0.04431094  0.58848592 -0.10070931

#bait = unlist(strsplit(c('SMC_EDNRA','SMC_FLT','EC_FLT1')[1],'\\_')) #,'EC_JCAD')
#infile = list.files('~/Desktop/MICOM/WhiteheadData/', full.names = T)[1]
#prepare(bait, infile)

#bait = unlist(strsplit(c('SMC_EDNRA','SMC_FLT','EC_FLT')[3],'\\_')) #,'EC_JCAD')
#infile = list.files('~/Desktop/MICOM/WhiteheadData/', full.names = T)[3]
#prepare(bait, infile)


