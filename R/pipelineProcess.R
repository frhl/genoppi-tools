#' @title pipeline for process proteomics data
#' @description a function that will process ip data
#' from different baits, files and cells
#' 

pipelineProcess <- function(baits, files, cells, steps = NULL){
  
  # standardize data
  baits = as.vector(baits); files = as.vector(files); cells = as.vector(cells)
  
  # 

}



baits = list('BLC2', 'MDM2', 'PTEN', 'TARDBP')
files = list()

