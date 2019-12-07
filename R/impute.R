#' @title Impute values
#' @description Replacing missing values with randomly sampled values from normal distribution,
#' with width SD x width and down-shifted Median-Sd x shift compared to observed sample distribution.
#' This is building upon the assumption that missing values have arisen due to low expression that 
#' can't be quantified. Therfore, shifting the median to lower expression levels will provide a 
#' proxy of this.
#' @param df a data.frame with numeric columns
#' @param stdwith numeric. change the factor of the standard deviation.
#' @param shift numeric. Negative values will shift the median distribution downwards.
#' @note No down-shifting and stdwith of 0.5 do not simualte low abudant missing values.
#' down-shifting of 0.8 and stdwidth of 0.5 simulates low abundant missing values. 
#' down-shifting of 3.6 and stdwith of 0.5 results in bi-modal distribution.
#' 
#' @references (Perseus, Tyanova et al. 2016)
#' @return data.frame with missing values imputed.


impute.gaussian <- function(df, stdwidth = 1, shift = 0){
  cols <- as.vector(unlist(lapply(df, function(x) is.numeric(x) & any(is.na(x)))))
  df$imputed <- as.numeric(apply(df, 1, function(x) any(is.na(x))))
  df[, cols] <- lapply(df[, cols], function(x){
    std <- sd(x, na.rm = T) * stdwidth # adjsuted/down-shifted mean (sample mean - SD * shift)
    med <- median(x, na.rm = T) + sd(x, na.rm = T) * shift  # adjsuted SD (sample SD * width)
    x[is.na(x)] <- rnorm(sum(as.numeric(is.na(x))), mean = med, sd = std)
    return(x)
  })
  return(df)
}




