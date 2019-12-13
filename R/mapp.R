#' @title map different gene names to one single gene name
#' @description Function for transferring a single gene or multiple genes to the consistent gene symbols
#' @param genes what genes should be in the query?
#' @param type protein or gene?
#' @param verbose print out prelimianry results such as, which were mapped? Did databases agree?
#' @note ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz was downloaded on 12DEC2019. This file was converted
#' into an inverted dict.
#' @examples \dontrun{mapp(LFS1)}
#' @export


mapp <- function(vec, type = 'gene', verbose = T){
 
  #load('ncbi_aliases')
  vec <- as.vector(vec)
  mapping = ncbi_aliases[[as.vector(vec)]]
  failed <- vec[is.na(mapping)]
  sucess <- vec[!is.na(mapping)]
  
  return(list(mapping = mapping, not.mapped = failed, was.mapped = sucess))  
   
}


#mapp('G3XAM7')
