# Install the released version from CRAN

.libPaths(c("~/Toolbox/rlib/",.libPaths()))
devtools::load_all()

test_that('normalize with median',{

  df <- data.frame(A=rnorm(100,5,15), B=runif(100,0,43))
  result = normalize(df, type = 'median')
  reference = as.data.frame(cbind(df$A-median(df$A), df$B-median(df$B)))
  colnames(reference) = c('A','B')
  expect_equal(result, reference)

})

test_that('normalize with mean',{
  
  df <- data.frame(A=rnorm(100,5,15), B=runif(100,0,43))
  result = normalize(df, type = 'mean')
  reference = as.data.frame(cbind(df$A-mean(df$A), df$B-mean(df$B)))
  colnames(reference) = c('A','B')
  expect_equal(result, reference)
  
})
