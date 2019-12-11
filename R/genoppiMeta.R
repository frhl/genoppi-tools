#' @title Genoppi databases
#' @param bait The bait that is to quries within the reference
#' @param reference the reference database. Currently only supports InWeb.
#' @family genoppi
#' @export

genoppi.database <- function(bait, reference = 'inweb'){
  
  require(hash)
  # generate a list of interactors based on the database
  # choosen. So far, only inweb is supported.
  if (tolower(reference) == 'inweb'){
    data(inweb_hash)
    inweb <- data.frame(gene=keys(inweb_hash))
    inweb$significant <- inweb$gene %in% inweb[[bait]]
    return(inweb)
  } else { stop(paste(reference, 'is not a valid database.')) }

}
  




