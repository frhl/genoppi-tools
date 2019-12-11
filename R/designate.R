#' @title Designate significant observations
#' @description Designate observations as significant based on filtering them using dplyr::filter.
#' @param data a data.frame
#' @param ... any arguments that can be inputted into dplyr::filter.
#' @author frhl
#' @examples \dontrun{
#' data = df %>% genoppi.ttest() ∂%>% designate(FDR < 0.1, pvalue < 0.05)}
#' @family genoppi
#' @export

designate <- function(data, ...){
  
  # check data.format
  cnames <- colnames(data)
  stopifnot('logFC' %in% cnames)
  stopifnot('FDR' %in% cnames)
  stopifnot('pvalue' %in% cnames)
  
  # Filter out the data needed and 
  # designate as significant
  data$idd <- 1:nrow(data)
  data$significant <- FALSE
  filtered = filter(data, ...)
  if (nrow(filtered) != 0) {
    data[data$idd %in% filtered$idd, ]$significant <- TRUE
  } else {warning('no rows were designated as significant since the criteria did not apply.') }
  data$idd <- NULL
  
  return(data)
  
}