
#' @title misc tools
#' @rdname misc
#' @export

#' @description length of unique
#' @export
lun <- function(x) length(unique(as.vector(x)))

#' @description returns true for x not in y
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
#' @export
warn <- function(msg){
  write(msg,stderr())
}
