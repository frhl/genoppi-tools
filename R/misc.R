#' @title misc tools
#' @rdname misc
#' @export

#' @description length of unique items.
#' @param x a vector or list of items.
#' @export
lun <- function(x) length(unique(as.vector(x)))

#' @description returns true for x not in y
#' @param x value x
#' @param y value y
#' @export
'%nin%' <- function(x, y) !(x %in% y)

#' @title omit nulls from list
#' @description remove NULLs in list
#' @export
null.omit <- function(x) {
  x[!vapply(x, is.null, logical(1))]
}

#' @title warnings to stderr
#' @description sends a message to stderr
#' @param msg the message.
#' @export
warn <- function(msg){
  write(msg,stderr())
}

#' @title make genoppi directories
#' @description Makes a directory according to bait, cell, imputation
#' @return will create a directory based on relative paths and return
#' a list of file names to be used for writing.
#' @export
mkdir <- function(bait, cell, impute){
  
  cur.time <- toupper(format(Sys.time(), "%d%b%Y"))
  newdir <- paste0(c('derived/',bait,'_', cell), collapse = '')
  output <- paste0(c(bait,'_', cell,'_IMPUTE=',impute,'_', cur.time), collapse = '')
  newfile <- file.path(newdir, output)
  dir.create(newdir)
  
  return(list(dir=newdir, output = output, basename = newfile))
}
