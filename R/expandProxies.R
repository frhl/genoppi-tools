#' @title Expand proxy interactions
#' @description Somtimes in Proteomics, the antibody is specific to another protein other
#' than the bait. This script will look throught the top enriched proteins (i.e. logFC >> 0),
#' and see whether they have interactors found in inweb and the current dataset.
#' @return a list of proteins that has inweb interactors in the dataset.
#' @export
#' @return a table that can be inputted to genoppi

expandProxies <- function(data, cutoff.fdr = 0.1, top=10, verbose = T){
  
  cnames = colnames(data)
  stopifnot('FDR' %in% cnames)
  stopifnot('logFC' %in% cnames)
  stopifnot('gene' %in% cnames)
  
  # Select the top most enriched data points
  tmp <- data[data$FDR < cutoff.fdr, ]
  tmp <- tmp[rev(order(tmp$logFC)), ]
  dat <- tmp[1:min(top, nrow(tmp)), ]
  
  # Find inweb interactors that are also present in data
  interact <- lapply(as.character(dat$gene), function(x){
    # subset of interactors already in the data
    proteins = interactors(x, verbose = F)
    proteins = proteins[proteins$significant,]
    proteins = proteins[proteins$gene %in% data$gene,]
    if (verbose) warn(paste('[InWeb]: Found', nrow(proteins), 'interactors for', x))
    return(proteins)
  })
  
  # return list of desired genes
  names(interact) <- as.character(dat$gene)
  return(interact)
  
}

#' @title Plot proxies sequentially
#' @description takes the data and a list of data.frames (generated with \code{expandProxies}),
#' and generated an overlap plot.
#' @param data a data.frame containing gene, logFC, pvalue and significant.
#' @param proxylist a list of data.frames in which the name of the list item correspond
#' to the prox bait, and the data.frame is all the interactors.

plotProxies <- function(data, proxylist, title.vs = 'control'){
  
  for (i in 1:length(proxylist)){
    bait = names(proxylist)[i]
    interactors = proxylist[[i]]
    title = paste(c('proxy bait:',bait, 'vs', title.vs, '(InWeb Overlap)'), collapse = ' ')
    plotOverlap(data, bait, interactors, title)
  }
}




