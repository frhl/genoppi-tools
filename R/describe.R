#' @title Describe data
#' @description Return boolean vector with length equal to number of columns, 
#' indicating which index contains intensity, ratios, baits column.
#' @param data a data.frame with proteomics data
#' @param verbose print to stdout
#' @author flassen
#' @export


describe <- function(data, control = 'mock', verbose = F){
  
  info <- list()
  # grab column details
  cnames <- tolower(colnames(data))
  info$cols.intensity<- grepl('int', cnames) # intensity columns
  info$cols.ratios <- grepl('ratio', cnames) # ratio columns
  info$cols.control <- grepl(tolower(control), cnames) & (!info$cols.ratios)
  info$possiblebaits <- info$cols.intensity & (!info$cols.ratios) & (!info$cols.control)
  
  # That accession column can be both accession
  info$col.accession <- grepl('accession', cnames)
  info$col.unique.proteins <- grepl('unique', cnames)
  
  return(info)
}
