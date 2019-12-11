#' @title Stamp a ggplot
#' @description Will stamp a ggplot in the lower right hand corner with the current time and date.
#' @param msg should there be a message before the date?
#' @param col color of the stamp.
#' @param alpha transparency
#' @param size size of the text
#' @export
#' @author Frederik

ggstamp <- function(msg = '', col = 'black', alpha = 0.9, size = 2){
  
  str = paste0(msg, format(Sys.time(), "%a %b %d %H:%M:%S %Y"))
  anno = data.frame(x=Inf,y=-Inf,hjust=1,vjust=0,label=str)
  p = geom_text(data=anno,aes(x=x, y=y, hjust=hjust, vjust=vjust, label=label), size = size, alpha = alpha, color = col)
  return(p)
}
