#' @title Genoppi databases of interactors
#' @description Input a bait and the current interactors found in the selected
#' database will return a list 
#' @param bait The bait that is to quries within the reference
#' @param reference the reference database. Currently only supports InWeb.
#' @param verbose how many interactors were found in the database?
#' @return a data.frame with genename and significant, that indicates whether
#' the bait-database interactor exists.
#' @family genoppi
#' @examples \dontrun{d1 <- interactors('CACNA1C')}
#' @export

interactors <- function(bait, verbose = F, reference = 'inweb'){
  
  require(hash)
  # generate a list of interactors based on the database
  # choosen. So far, only inweb is supported.
  if (tolower(reference) == 'inweb'){
    data(inweb_hash)
    inweb <- data.frame(gene=keys(inweb_hash))
    inweb$significant <- inweb$gene %in% inweb_hash[[bait]]
    if (verbose) warn(paste('[InWeb]: Found', sum(inweb$significant), 'interactors for', bait))
    return(inweb)
  } else { stop(paste(reference, 'is not a valid database.')) }

}






