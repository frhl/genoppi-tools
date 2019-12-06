#' @title Describe data
#' @description A function that will take a data.frame
#' and describe the columns, type of data which will be returned in
#' a list
#' @param data a data.frame with proteomics data
#' @param verbose print to stdout
#' @author flassen
#' @export


describe <- function(data, verbose = F){
  
  info <- list()
  cnames <- tolower(colnames(data))
  info$cols.intensity<- grepl('int', cnames)
  info$cols.ratios <- grepl('ratio', cnames)
  info$cols.control <- grepl('mock', cnames) & (!info$cols.ratios)
  info$possiblebaits <- info$cols.intensity & (!info$cols.ratios) & (!info$cols.control)
  info$col.accession <- grepl('accession', cnames)
  info$col.unique.proteins <- grepl('unique', cnames)
  
  # check format of columns
  if (verbose){
    if (!all(sapply(data[, (info$cols.intensity)], is.numeric))) write('note: intensity columns are not numeric.',stderr())
    if (!all(sapply(data[, (info$cols.control)], is.numeric))) write('note: mock columns are not numeric.',stderr())
  }
  # check experiment
  if (grepl('itraq', paste(cnames, collapse = ''))) info$experiment <- 'iTRAQ'
  return(info)
}
