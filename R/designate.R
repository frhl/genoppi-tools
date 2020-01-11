#' @title Designate significant observations
#' @description Designate observations as significant based on filtering them using dplyr::filter.
#' Can be used in conjunction with genoppi plotting functions.
#' @param data a data.frame
#' @param ... any arguments that can be inputted into dplyr::filter.
#' @author frhl
#' @examples \dontrun{
#' data = df %>% genoppi.ttest() ∂%>% designate(FDR < 0.1, pvalue < 0.05)}
#' @family genoppi
#' @export

designate <- function(data, ...){
  
  require(dplyr)
  
  cnames <- colnames(data)
  stopifnot('logFC' %in% cnames)
  stopifnot('FDR' %in% cnames)
  stopifnot('pvalue' %in% cnames)
  data$idd <- 1:nrow(data)
  data$significant <- FALSE
  filtered = filter(data, ...)
  if (nrow(filtered) != 0) {
    data[data$idd %in% filtered$idd, ]$significant <- TRUE
  }
  data$idd <- NULL
  
  return(data)
  
}