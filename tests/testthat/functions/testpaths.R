#' @title Get resulting and referencing test files
#' @param id id of the test, e.g. 'A1'
#' @param fun name of the function/test as per file saved.


testpaths <- function(id, fun){
  core = 'tests/testthat'
  alt = c('reference', 'results')
  paths = lapply(alt, function(x){
    path = paste0(c(core, x, fun, ''), collapse = '/')
    files = list.files(path, pattern = id, full.names = T)
    if (length(files) == 0) files <- NULL
    return(files)
  })
  names(paths) <- c('ref','res')
  return(paths)
}
