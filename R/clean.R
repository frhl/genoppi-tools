
#' @title clean
#' @description filter out contaminants and proteins with
#' less than < 2 unique peptides


clean <- function(table, protein_column = NULL, threshold = 0, verbose = TRUE){
  
  # assume first column if none is specified
  if (is.null(protein_column)) protein_column = colnames(table)[1]
  if (!is.character(table[, protein_column])) stop(paste(protein_column, 'does not contain protein/gene names'))
  
  # get numeric columns
  columns_numeric <- sapply(table, is.numeric)
  table_numeric <- table[, columns_numeric]
  
  # search rowwise
  #rows_ok <-  as.logical(apply(sapply(table_numeric, function(x) x < 2), 1, sum))
  #write(paste(length(rows_ok) - sum(rows_ok), 'rows was removed.'),stderr())
  #table <- table[rows_ok, ]
  
  
  return(table)
}
