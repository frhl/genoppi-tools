#' @title Calculate a moderated t-test
#' @description alculate a moderated t-test
#' @param df a data.frame with 'rep' in column names.
#' @author April Kim
#' @export


moderatedTTest <- function(df){
  
  calculated <- df
  colnames(calculated)[1] <- "gene"
  
  # apply median normalization
  medianNorm <- function(d){return(d - median(d, na.rm = TRUE))}
  for (i in 2:dim(df)[2]) {
    colStr <- paste("rep",i-1,sep="")
    calculated[,colStr] <- medianNorm(calculated[,colStr])
  }
  
  # calc moderated t-test
  myfit <- lmFit(calculated[,-1], method="robust")
  myfit <- eBayes(myfit)
  modtest <- topTable(myfit, number=nrow(myfit), sort.by='none')
  colnames(modtest)[4:5] <- c("pvalue","FDR")
  
  # return data frame with test results: gene, rep1, rep2, ..., logFC, pvalue, FDR 
  calculated <- data.frame(cbind(calculated, modtest[,-c(2,3,6)]))
  calculated <- calculated[with(calculated, order(-logFC, FDR)),]
  
  return(calculated)
}