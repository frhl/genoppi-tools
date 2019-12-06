#' @title Transform data
#' @description Normalize a table of label or label free data. 
#' The function will transform all numeric columns.
#' intensity value.
#' @param table a data.frame containing numeric columns
#' @param type an r-command like 'log2' or 'log10'.
#' @author flassen
#' @export

pTransform <- function(table = NULL, type = 'log2'){
  
  columns_numeric <- sapply(table, is.numeric)
  table_transformed <- do.call(type, list(table[, columns_numeric]))
  table[, columns_numeric] <- table_transformed
  return(table)
  
}

