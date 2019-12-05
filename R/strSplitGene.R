#' @title splits an accession id
#' @description splits the gene in a column/vector of accession id
#' so that the gene is extract.
#' @author flassen
#' @export

strSplitGene <- function(vec){
  # modify accession to get gene name
  gene = unlist(lapply(strsplit(as.character(vec), '\\|'), function(x) x[length(x)]))
  gene = unlist(lapply(strsplit(gene, '\\_'), function(x) x[1]))
  if (any(is.na(gene))) warn('warning. some accession numbers correctly processed')
  return(gene)
}