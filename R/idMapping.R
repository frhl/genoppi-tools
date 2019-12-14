#' @title splits an accession id
#' @description splits the gene in a column/vector of accession id
#' so that the gene is extract.
#' @author flassen
#' @export

strSplitGene <- function(vec, get = 'symbol'){
  # modify accession to get gene name
  browser()
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
    entry <- unlist(strsplit(as.character(x),'\\||\\_'))
    if (length(entry) == 3) return(c(NA, entry))
    else if (length(entry == 4)) return(entry)
    else return(rep(NA,4))
  })
  
  
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
  hits = hm$keys() %in% mat$uniprot$id
  
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






#uni <- as.character(mat$uniprot$id) %in% dat$uniprot_gn_id
#sym <- as.character(mat$symbol) %in% dat$hgnc_symbol

#ncbi_aliases[['U3KPZ7']]


#uniprotid_to_hgnc <- hashmap(keys=hgnc_uniprot_mapping$uniprot_gn_id, values=hgnc_uniprot_mapping$hgnc_symbol)
#save_hashmap(uniprotid_to_hgnc, file = 'uniprotid_to_hgnc')


#.libPaths('~/Toolbox/rlib/')
#library('biomaRt')
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
#dat = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'uniprot_gn_symbol'), mart = ensembl)
dat = getBM(attributes = c('hgnc_symbol', 'uniprot_gn_id', 'uniprot_gn_symbol'), mart = ensembl)

#hgnc_uniprot_mapping = dat
#save(hgnc_uniprot_mapping, file = 'hgnc_uniprot_mapping.RData')

#htable <- mkhashtable(dat)

