#' @title Expand proxy interactions
#' @description Somtimes in Proteomics, the antibody is specific to another protein other
#' than the bait. This script will look throught the top enriched proteins (i.e. logFC >> 0),
#' and see whether they have interactors found in inweb and the current dataset.
#' @return a list of proteins that has inweb interactors in the dataset.
#' @export
#' @return a table that can be inputted to genoppi

expandProxies <- function(data, cutoff.fdr = 0.1, top=2, verbose = T){
  
  # Select the top most enriched data points
  tmp <- data[data$FDR < cutoff.fdr, ]
  tmp <- tmp[rev(order(tmp$logFC)), ]
  dat <- tmp[1:top, ]
  
  # Find inweb interactors that are also present in data
  interact <- lapply(as.character(dat$gene), function(x){
    proteins = interactors(x, verbose = F)
    proteins = proteins[proteins$significant,]
    proteins = proteins[proteins$gene %in% data$gene,]
    if (verbose) warn(paste('[InWeb]: Found', nrow(proteins), 'interactors for', x))
    return(proteins)
  })
  
  # return list of desired results
  names(interact) <- as.character(dat$gene)
  return(interact)
  
}


