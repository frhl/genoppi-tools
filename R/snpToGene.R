#' @title map snps to a specific gene
#' @description maps a snp to a gene by using a dict as reference
#' @param snp a vector of snps 
#' @author flassen
#' @export
#' @note uses \code{null.omit} from misc.R
#' @examples \dontrun{snpToGene('rs142920272')}


snpToGene <- function(snp){

  if ('genes_snps' %nin% ls(envir = .GlobalEnv)) 
    {write('reading SNP/Gene refrence..',stderr()); data('genes_snps')}
  genes = names(genes_snps)
  result = lapply(genes, function(g){
    table_snps = genes_snps[[g]]
    check = table_snps %in% as.vector(snp) 
    if (any(check)){return(table_snps[check])}
  })
  names(result) = genes
  return(null.omit(result))
}





