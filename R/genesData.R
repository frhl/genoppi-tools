#' gene table (old)
#'
#' data from genoppi used to map snps to gene. 
#' Todo: ask april about reference, think its hg19
#'
#' @docType data
#'
#' @usage data(genes_snps)
#'
#'
#' @keywords datasets
#'
#'
#'
#' @examples
#' data(genes_snps)
#' 
"genes_snps"


#' InWeb Hash Table
#'
#' InWeb Hash table downloaded on october 2018.
#' 
#' @docType data
#' @author Yu-Han Shu
#'
"inweb_hash"


#' NCBI Hash Table for mapping genes
#'
#' @description  NCBI inverted dictionary hash table that maps an alias to the same gene format.
#' @note ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz was downloaded on 12DEC2019. This file was converted
#' into an inverted dict. Was downloaded according to this: https://www.biostars.org/p/1378/
#' @docType data
#' @author Frederik Lassen
#'
"ncbi_aliases"


#' Uniprot for mapping protein names to primary name
#'
#' @description  Uniprot inverted dictionary hash table that maps an alias to the same protein format.
#' @note Downloaded 12DEC2019 by https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+proteome%3Aup000005640,
#' and selected for 'gene names (synonyms) and 'gene names (primary)'. The synonyms provides the keys and the primary
#' provide the corresponding value.
#' @docType data
#' @author Frederik Lassen
#'
"uniprot_aliases"

#' map all entries to ensemble ID
#' 
#' @description Maps all the indicates entries to ensemble ID.
#' @note does
#' 
#' 


