#' @title Calculate a moderated t-test
#' @description alculate a moderated t-test
#' @param df a data.frame with 'rep' in column names.
#' @author April Kim
#' @export


moderatedTTest <- function(df){
  
  
  # isolate rep columns
  reps <- grepl('rep', colnames(df)) & unlist(lapply(df, is.numeric))
  calculated <- df[, reps]
  
  # apply median normalization
  # note: is this step needed, since intensity values
  # have already been normalized?
  calculated <- normalize(calculated)
  
  # calc moderated t-test
  myfit <- lmFit(calculated, method="robust")
  myfit <- eBayes(myfit)
  modtest <- topTable(myfit, number=nrow(myfit), sort.by='none')
  colnames(modtest)[4:5] <- c("pvalue","FDR")
  
  # return data frame with test results: gene, rep1, rep2, ..., logFC, pvalue, FDR 
  calculated <- data.frame(cbind(calculated, modtest[,-c(2,3,6)]))
  calculated <- calculated[with(calculated, order(-logFC, FDR)),]
  calculated <- cbind(df$gene, calculated)
  colnames(calculated)[1] <- 'gene'
  
  return(calculated)
}