#' @title Calculate a moderated t-test
#' @description Will do a moderated t.test pairwise on replicate columns, and
#' calculate the FDR and p-value.
#' @param df a data.frame with 'rep' in column names.
#' @param keep a vector of characters containing colnames that should be kept in the data.frame
#' @author April Kim / Frederik Heymann
#' @family genoppi
#' @export
 

mttest <- function(df, keep = c('imputed')){
  
  require(limma)
  
  # check input
  stopifnot('gene' %in% colnames(df))
  
  # isolate rep columns
  reps <- grepl('rep', colnames(df)) & unlist(lapply(df, is.numeric))
  calculated <- df[, reps]

  
  # calc moderated t-test
  myfit <- lmFit(calculated, method="robust")
  myfit <- eBayes(myfit)
  modtest <- topTable(myfit, number=nrow(myfit), sort.by='none')
  colnames(modtest)[4:5] <- c("pvalue","FDR")
  
  # return data frame with test results: gene, rep1, rep2, ..., logFC, pvalue, FDR 
  calculated <- data.frame(cbind(calculated, modtest[,-c(2,3,6)]))
  calculated <- cbind(df$gene, calculated)
  colnames(calculated)[1] <- 'gene'
  
  # keep columns
  keep = keep[(keep %in% colnames(df))]
  ndf = data.frame(df[,keep])
  names(ndf) = keep
  calculated <- cbind(calculated, ndf)
  
  # order columns
  calculated <- calculated[with(calculated, order(-logFC, FDR)),]

  return(calculated)
}