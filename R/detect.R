#' @title Detect protein and gene IDs
#' @description Will look up acession IDs to
#' see whether the entry is in it.
#' @param data a data.frame
#' @param entry a string or regular expression.
#' @author flassen
#' @family id
#' @return row-wise boolean indiciating whether the entry is present in the data.
#' @note
#' @export


detect <- function(data, entry){
  
  if (!is.null(entry)){
    if (is.character(data)) data <- read.csv(data)
    if (is.vector(entry)) entry <- paste(entry,collapse='|')
    cols <- describe(data)
    mata <- acession.matrix(data[,cols$col.accession])
    matc <- acession.convert(mata, verbose = F)
    data <- cbind(data, matc)
    ## get row in whitch the entry can be found
    result <- apply(sapply(matc, function(x) grepl(entry, x)), 1, any) 
    return(result)
  } else return(FALSE)

}
