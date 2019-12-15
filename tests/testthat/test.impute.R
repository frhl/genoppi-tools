

test_that('imputation simulating low abundant proteins',{

  ### n is big
  set.seed(1337)  
  df <- data.frame(A=rnorm(10000, 0, sd = 4), B = rnorm(10000, 1, sd = 3))
  sam <- sample(1:10000, 1000)
  df[sam, 1] <- NA
  df <- impute.gaussian(df, width = 0.3, shift = -1.8)
  x1 <- hist(df[df$imputed==1, ]$A, 100, xlim = c(-15, 15))
  x2 <- hist(df[df$imputed==0, ]$A, 100, xlim = c(-15, 15))
  
  # visual inspection
  #plot(x2, col = 'grey', main = 'desireable imputation') # Plot 1st histogram using a transparent color
  #plot(x1, col = 'yellow', add = TRUE) # Add 2nd histogram using different color
  
  # expect lower mean
  expect_equal(mean(df[df$imputed==1, ]$A), -7.318569, tolerance = 0.001)
  expect_equal(mean(df[df$imputed==0, ]$A), -0.02240977, tolerance = 0.001)

  ## n is small
  set.seed(1337)  
  df <- data.frame(A=rnorm(2000, 0, sd = 4), B = rnorm(2000, 1, sd = 3))
  sam <- sample(1:2000, 100)
  df[sam, 1] <- NA
  #df[sam, 1] <- NA
  df <- impute.gaussian(df, width = 0.5, shift = -1.8)
  x1 <- hist(df[df$imputed==1, ]$A, 50, xlim = c(-15, 15))
  x2 <- hist(df[df$imputed==0, ]$A, 50, xlim = c(-15, 15))
  
  # visual inspection
  plot(x2, col = 'grey', main = 'desireable imputation') # Plot 1st histogram using a transparent color
  plot(x1, col = 'yellow', add = TRUE) # Add 2nd histogram using different color
  
  
  
})


test_that('imputation without shifting',{
  
  set.seed(1337)  
  df <- data.frame(A=rnorm(10000, 0, sd = 4), B = rnorm(10000, 1, sd = 3))
  sam <- sample(1:10000, 1000)
  df[sam, 1] <- NA
  df <- impute.gaussian(df, width = 0.5, shift = 0)
  x1 <- hist(df[df$imputed==1, ]$A, 100, xlim = c(-15, 15))
  x2 <- hist(df[df$imputed==0, ]$A, 100, xlim = c(-15, 15))
  
  # visual inspection
  plot(x2, col = 'grey') # Plot 1st histogram using a transparent color
  plot(x1, col = 'yellow', add = TRUE) # Add 2nd histogram using different color
  
  # expect approx same mean
  expect_equal(mean(df[df$imputed==1, ]$A), -0.08632506, tolerance = 0.0001)
  expect_equal(mean(df[df$imputed==0, ]$A), -0.02240977, tolerance = 0.0001)
  
})







