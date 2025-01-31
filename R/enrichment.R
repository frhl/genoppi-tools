#' @title Calculate enrichment using Fisher exact test
#' @description Calculate enrichment of a data and reference data base for
#' a specific bait and its interactors. Will do p-test using a fisher exact test
#' based on the overlap between the main data and the reference data which both
#' interacts with bait.
#' @param data data.frame containing columns gene and significance
#' @param bait designated bait in the data.
#' @param reference data base reference
#' @family genoppi
#' @export


enrichment <- function(data, bait, reference){
  
  ## check format
  cnames = colnames(data)
  stopifnot('significant' %in% colnames(data))
  stopifnot('gene' %in% colnames(data))
  stopifnot('gene' %in% colnames(reference))
  
  ## ready data for fisher exact test
  preys <- unique(data$gene[data$gene != bait])
  
  # significant preys
  sigPreys <- unique(data$gene[data$significant & data$gene %in% preys])		
  
  # significant genes in gene list, limited to detected preys
  sigGenes <- reference$gene[reference$significant & reference$gene %in% preys]
  
  # overlap of significant preys and significant genes
  overlap <- intersect(sigPreys,sigGenes)
  # sig preys but not sig genes
  sigPreysOnly <- setdiff(sigPreys,sigGenes)
  # sig genes but not sig preys
  sigGenesOnly <- setdiff(sigGenes,sigPreys)
  # neither sig preys nor sig genes
  neither <- setdiff(setdiff(preys,sigPreys),sigGenes)
  
  # Fisher's exact test (one-tailed)
  fisherM <- matrix(c(length(overlap),length(sigPreysOnly),
                      length(sigGenesOnly),length(neither)),nrow=2)
  fisherP <- fisher.test(fisherM,alternative="greater")$p
  
  return(list(sigGenes = sigGenes, sigPreys = sigPreys, overlap = overlap, fisherP = fisherP))
}


enrichment_inweb <- function(data, bait, reference){
  
  ## check format
  cnames = colnames(data)
  stopifnot('significant' %in% colnames(data))
  stopifnot('gene' %in% colnames(data))
  stopifnot('gene' %in% colnames(reference))
  
  ## ready data for fisher exact test
  preys <- unique(data$gene[data$gene != bait & data$gene %in% reference$gene])
  
  # significant preys
  sigPreys <- unique(data$gene[data$significant & data$gene %in% preys])		
  
  # significant genes in gene list, limited to detected preys
  sigGenes <- reference$gene[reference$significant & reference$gene %in% preys]
  
  # overlap of significant preys and significant genes
  overlap <- intersect(sigPreys,sigGenes)
  # sig preys but not sig genes
  sigPreysOnly <- setdiff(sigPreys,sigGenes)
  # sig genes but not sig preys
  sigGenesOnly <- setdiff(sigGenes,sigPreys)
  # neither sig preys nor sig genes
  neither <- setdiff(setdiff(preys,sigPreys),sigGenes)
  
  # Fisher's exact test (one-tailed)
  fisherM <- matrix(c(length(overlap),length(sigPreysOnly),
                      length(sigGenesOnly),length(neither)),nrow=2)
  fisherP <- fisher.test(fisherM,alternative="greater")$p
  
  return(list(sigGenes = sigGenes, sigPreys = sigPreys, overlap = overlap, fisherP = fisherP))
}


fisher.pvalue.inweb <- function(bait, cell_data_all, cell_data_significant){
  
  preys <- unique(intersect(hash::keys(inweb_hash),cell_data_all))  #unique(inweb_data[inweb_data %in% cell_data_all & inweb_data %nin% bait])
  preys <- preys[preys %nin% bait]
  # significant preys
  #sigPreys <- unique(data$gene[data$significant & data$gene %in% preys])	
  sigPreys <- cell_data_all[cell_data_all %in% cell_data_significant & cell_data_all %in% preys]
  # significant genes in gene list, limited to detected preys
  #sigGenes <- reference$gene[reference$significant & reference$gene %in% preys]
  sigGenes <- inweb_data[inweb_data %in% preys]
  # overlap of significant preys and significant genes
  overlap <- intersect(sigPreys,sigGenes)
  # sig preys but not sig genes
  sigPreysOnly <- setdiff(sigPreys,sigGenes)
  # sig genes but not sig preys
  sigGenesOnly <- setdiff(sigGenes,sigPreys)
  # neither sig preys nor sig genes
  neither <- setdiff(setdiff(preys,sigPreys),sigGenes)
  
  fisherM <- matrix(c(length(overlap),length(sigPreysOnly),
                      length(sigGenesOnly),length(neither)),nrow=2)
  fisherP <- fisher.test(fisherM,alternative="greater")$p
  
  return(fisherP)
  
}


