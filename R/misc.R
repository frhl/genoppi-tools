#' @rdname misc
#' @export

#' @description length of unique
#' @export
lun <- function(x) length(unique(as.vector(x)))

#' @description returns true for x not in y
#' @export
'%nin%' <- function(x, y) !(x %in% y)

#' @description remove NULLs in list
#' @export
null.omit <- function(x) {
  x[!vapply(x, is.null, logical(1))]
}
