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
    # subset of proteins already in the data
    proteins = interactors(x, verbose = F)
    proteins = proteins[proteins$significant,]
    proteins = proteins[proteins$gene %in% data$gene,]
    #proteins = data[data$gene %in% proteins$gene,]
    #proteinsSorted <- proteins[rev(order(proteins$logFC)), ]
    #proteinsSorted$significant <- TRUE # all are found in inweb
    
    if (verbose) warn(paste('[InWeb]: Found', nrow(proteins), 'interactors for', x))
    return(proteins)
  })
  
  # return list of desired results
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

makeInteractionMap <- function(data, ...){
  
  stop('work in progress..')
  
  dat <- filter(data, ...)
  mapping <- lapply(as.character(dat$gene), function(x){
    proteins = interactors(x, verbose = T)
    return(proteins)
  })
  
  # data into inweb matrix
  genes <- keys(inweb_hash)
  comb <- lapply(mapping, function(x) x[2])
  comb <- do.call(cbind, comb)
  colnames(comb) <-  as.character(dat$gene)
  rownames(comb) <- comb$genes
  
  # filter data
  indata <- comb[comb$genes %in% dat$gene, ] # note: here we take all the proteins
  indata[, ]
  
  #indata[, 2:ncol(indata)] <- indata[, c('genes', as.character(indata$genes))]
  
}





