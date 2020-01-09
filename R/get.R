

## remove these functions at some point

get.interactors.genes <- function(data, known.interactors){
  
  known.interactors.lst <- as.character(known.interactors[known.interactors$significant == TRUE, ]$gene)
  in.data = data[data$gene %in% known.interactors.lst, ]
  interactors.genes <- NA
  if (nrow(in.data) > 0){
    interactors.ok <- designate(in.data, FDR < 0.1, logFC > 0)
    if (any(interactors.ok$significant)){
      interactors.genes <- paste(interactors.ok[interactors.ok$significant == TRUE, ]$gene, collapse=';')
    }
  } 
  return(interactors.genes)

}


## get all interactars in the dataset, regardless of significane
get.indata.interactors <- function(data, known.interactors){
  
  known.interactors.lst <- as.character(known.interactors[known.interactors$significant == TRUE, ]$gene)
  in.data = data[data$gene %in% known.interactors.lst, ]
  interactors.genes <- paste(in.data$gene, collapse=';')
  return(interactors.genes)

}



get.summary.aggregate <- function(directory){
  
  file = list.files(directory, recursive = T, full.names = T, pattern = '8JAN.+SUMMARY')
  tab <- do.call(rbind, lapply(file, function(x) fread(x)))
  
  
}