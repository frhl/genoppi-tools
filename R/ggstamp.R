#' @title Stamp a ggplot
#' 


ggstamp <- function(msg = '', col = 'black', alpha = 0.9, size = 2){
  
  str = paste0(msg, format(Sys.time(), "%a %b %d %H:%M:%S %Y"))
  anno = data.frame(x=Inf,y=-Inf,hjust=1,vjust=0,label=str)
  p = geom_text(data=anno,aes(x=x, y=y, hjust=hjust, vjust=vjust, label=label), size = size, alpha = alpha, color = col)
  return(p)
}
