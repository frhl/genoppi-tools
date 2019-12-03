#' @title table genes and snps
#' @description creates a 2xn table with genes in the first
#' column and the associated snps in the second column.
#' @author flassen
#' @export
#' @examples \dontrun{geneTable(snp2gene('rs142920272'))}

geneTable <- function(snpgene, outfile = NULL, delim = ';'){

  dat = lapply(snpgene, function(x) paste(x, collapse = delim))
  df = data.frame(geneName = names(dat), snps = unlist(dat))
  row.names(df) <- NULL
  if (!is.null(outfile)) write.csv(df, outfile)
  return(df)
}

