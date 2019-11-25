#' @title Normalize data
#' @description Normalize a table of label or label free data. 
#' The function will calculate the \code{type} e.g. median of
#' all numeric columns and subtract the median from the 
#' intensity value.
#' @author flassen
#' @export

# for testthat: table <- data.frame(A=rnorm(100,5,15), B=runif(100,0,43), C = rep(c('A',"B"),50))


normalize <- function(table = NULL, type = 'median'){

  columns_numeric <- sapply(table, is.numeric)
  table_numeric <- table[, columns_numeric]
  values_type <- sapply(table_numeric, function(x) do.call(type, list(x)))
  table_normalized <- table_numeric - values_type
  table[, columns_numeric] <- table_normalized
  return(table)
  
}












