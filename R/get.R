

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
  
  #

  file = list.files('~/Projects/03_MICOM/derived/', recursive = T, full.names = T, pattern = '9JAN.+SUMMARY')
  
  library(data.table)
  
  micom_summary <- lapply(unique(file), function(x) read.csv(unique(x))[1,])
  micom_summary
  
  tab <- do.call(rbind, micom_summary)
  
  
  
  tab$martin <- gsub('03_MICOM','',tab$file)
  tab$martin <- as.numeric(gsub('martin','',regmatches(tab$martin,regexpr("[0-9]+martin",tab$martin))))
  tab <- tab[order(tab$martin), ]
  boxplot(tab$rep.cor~tab$martin, main = 'variability of replicate correlations')
  tab$data.bait.enriched[is.na(tab$data.bait.enriched)] <- FALSE

  write.csv(tab, file = '~/Projects/03_MICOM/analysis.summary.csv')
  
  #tab$data.bait.enriched <- as.character(tab$data.bait.enriched)
  #tab$data.bait.enriched[tab$data.bait.enriched == 'FALSE'] <- 'BLACK'
  #tab$data.bait.enriched[tab$data.bait.enriched == 'TRUE'] <- 'RED'
  #plot(tab$rep.cor~tab$martin, col = tab$data.bait.enriched)
  
  sum(tab$data.bait.enriched)/41
  
  
}