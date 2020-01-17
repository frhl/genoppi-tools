#' @title splits an accession id
#' @description splits the gene in a column/vector of accession id
#' so that the gene is extract.
#' @author flassen
#' @export

strSplitGene <- function(vec, get = 'symbol'){
  # modify accession to get gene name
  warn('deprecated. Will be reomved in future release.')
  if (get %in% 'symbol'){
    gene = unlist(lapply(strsplit(as.character(vec), '\\|'), function(x) x[length(x)]))
    gene = unlist(lapply(strsplit(gene, '\\_'), function(x) x[1]))
    if (any(is.na(gene))) warn('warning. some accession numbers correctly processed')
  } else if (get %in% 'id'){
    
    
  }
  return(gene)
}

#' @title Split acession ID into matrix
#' @description Generase an acession ID matrix that
#' can be further subsetted according to isoforms, gene symbol
#' and species. 
#' @author flassen
#' @family id
#' @export

acession.matrix <- function(vec){
  
  ## inital checks
  stopifnot(!is.null(vec))
  stopifnot(length(vec) > 0)

  ## assumes format X|X|X_X, i.e.  4 elements
  mat <- lapply(vec, function(x) {
    
    #if (grepl('FLAG', x)) browser()
    entry = unlist(strsplit(as.character(x),'\\||\\_'))
    n = as.numeric(length(entry))
    #if (grepl('FLAG', x)) {warning('FLAG ')} 
    if (n == 3) {return(c(NA, entry))}
    if (n == 4) {return(entry)}
    if (n == 5) {return(c(NA, NA, entry[4], NA))} # gi|999627|pdb|1EPT|B
    if (n == 6) {return(entry[3:6])} # for entries looking like this: gi|125138|sp|P01840.1|KAC4_RABIT
    return(rep(NA,4))
  })
  
  #unique(sort(unlist(lapply(mat, length))))
  #num <- which(lapply(mat, length) %in% c(5,6))
  #vec[num]
  
  mat <- as.data.frame(do.call(rbind, mat))
  colnames(mat) <- c('sp', 'uniprot', 'symbol', 'species')
  
  ## deal with isoforms
  isoforms <- lapply(strsplit(as.character(mat$uniprot), '\\-'), function(x){
    if (length(x) == 1) return(c(x, NA)) else return(x)})
  isoforms <- as.data.frame(do.call(rbind, isoforms))
  colnames(isoforms) <- c('uniprot.id', 'uniprot.isoform')
  
  ## check formats
  stopifnot(nrow(isoforms) == nrow(mat))
  
  ## merge data.frames
  mat <- as.data.frame(cbind(mat, isoforms))
  mat$uniprot <- NULL # remove old uniprot
  
  return(mat)
}

#' @title Map uniprot ID to HGNC symbol
#' @description Input an acession matrix to subset it accordingly. It will
#' map uniprot IDs to their corresponding HGNC symbols.
#' @param mat an acession.matrix
#' @author flassen
#' @family id
#' @return an acession matrix with two new columns:
#' 1) from will indicate whether the uniprot ID was found and converted to HGNC, NA will indicate this
#' was not the case (i.e. the original symbol found in the acession string will be used).
#' 2) hgnc the symbol derived after conversion from uniprot ID to HGNC.
#' @export

acession.convert <- function(mat,  verbose = T){
  
  require(hashmap)

  # check input
  stopifnot(!is.null(mat))
  stopifnot(length(dim(mat)) > 1)
  hm <- load_hashmap('~/Toolbox/packages/pRoteomics/data/uniprotid_to_hgnc')
  hits = hm$keys() %in% mat$uniprot.id
  
  ## assign new IDs
  mat$from <- NA
  mat$hgnc <- NA
  newids <- hm[[mat$uniprot.id]]
  mat$hgnc[!is.na(newids)] <- newids[!is.na(newids)]
  mat$hgnc[mat$hgnc == ''] <- NA
  mat$from[!is.na(mat$hgnc)] <- 'uniprot'
  mat$hgnc[is.na(mat$from)] <- as.character(mat$symbol[is.na(mat$from)])
  
  ## output
  coverage <- round((sum(hits)/nrow(mat))*100,2)
  if (verbose) warn(paste0('[uniprot] ', sum(hits),'/',nrow(mat),' (',coverage,'%) coverage of acession IDs.' ))
  return(mat)
  
}
