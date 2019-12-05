#' @title Transform data
#' @description Calculates log fold change by substracting
#' a case column with a bait column
#' @author flassen
#' @export

logFC <- function(df){
  # assummes log2(bait1), log2(control1), log2(bait2).. 
  pairs <- (ncol(df)-1)/2
  dfNum <- sapply(df, is.numeric)
  dfNum <- df[, dfNum]
  for (i in 1:pairs){
    colMock = dfNum[, (i*2)]
    colBait = dfNum[, (i*2)-1]  # log2(bait) - log2(control)
    dfNum[[paste0('rep',i)]] = colBait - colMock
  }
  return(cbind(df, dfNum[,grepl('rep', colnames(dfNum))]))
}