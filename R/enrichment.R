#' @title Calculate enrichment using Fisher exact test
#' @description 
#' @param data data.frame containing columns gene and significance
#' @param bait designated bait in the data.
#' @param reference data base reference
#' @family genoppi
#' 


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
