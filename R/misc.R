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
#' @param bait what bait was used?
#' @param cell what cell type?
#' @param name the name that is to be saved. Default is a combination of bait and cell.
#' @param variation add a variation to the outputted name
#' @return will create a directory based on relative paths and return
#' a list of file names to be used for writing.
#' @export
mkdir <- function(bait='', cell='', variation = "", run = '', name = paste0(bait,'_',cell)){
  
  # if run is specified, it will create a master folder
  # in which all files will be organized
  if (run != ''){
    dir.run = paste0('derived/',run, '/')
    dir.create(dir.run)
    newdir <- paste0(c(dir.run,name), collapse = '/')
  } else{
    newdir <- paste0(c('derived/',name), collapse = '')
  }
  cur.time <- toupper(format(Sys.time(), "%d%b%Y"))
  output <- paste0(c(name,'_',variation, cur.time), collapse = '')
  newfile <- file.path(newdir, output)
  pdfpath <- paste0(newfile, '.pdf')
  txtpath <- paste0(newfile, '.txt')
  csvpath <- paste0(newfile, '.csv')
  dir.create(newdir)
  
  return(list(dir=newdir, output = output, basename = newfile, pdfpath = pdfpath, txtpath = txtpath, csvpath = csvpath))
}
