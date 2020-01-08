#' @title Describe data
#' @description A function that will take a data.frame
#' and describe the columns, type of data which will be returned in
#' a list
#' @param data a data.frame with proteomics data
#' @param verbose print to stdout
#' @author flassen
#' @export


describe <- function(data, control = 'mock', verbose = F){
  
  info <- list()
  cnames <- tolower(colnames(data))
  info$cols.intensity<- grepl('int', cnames)
  info$cols.ratios <- grepl('ratio', cnames)
  info$cols.control <- grepl(tolower(control), cnames) & (!info$cols.ratios)
  info$possiblebaits <- info$cols.intensity & (!info$cols.ratios) & (!info$cols.control)
  
  # That accession column can be both accession
  # or gene id
  #if (any(grepl('accession', cnames))){
    info$col.accession <- grepl('accession', cnames)
  #} else info$col.accession <- grepl('^gene', cnames)

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
