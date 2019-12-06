#' @title pipeline for process proteomics data
#' @description a function that will process ip data
#' from different baits, files and cells
#' @param bait a vector containing that contains the bait and cell type that should be matched in a column.
#' @param infile the file that contains the raw data, i.e. accession numbers, intensity values, ratios etc.
#' @param cols optional manual entry. A vector of columns that are present in the dataset header. Follows 
#' the format of columns: acession, bait1, mock1, bait2, mock2, bait3, ..
#' @param impute how should missing data be imputed? NULL means that missing data rows are dropped.
#' @param transform charcacter. an R-command for how the data should be transformed.
#' @param normalization character. an R-command for how the data should be collumn-wise transformned.
#' @param filter character. only accession IDs of the filter specified are included.
#' @param verbose boolean. if true, returns the table and a list with updates.
#' @return a table that can be inputted to genoppi

prepare <- function(bait, infile, cols = NULL, impute = NULL, transform = 'log2', normalization = 'median', filter = "HUMAN", verbose = F){
  
  # check input
  if (any(!file.exists(infile))) stop('one of the inputted files does not exist.')

  ## do some initial checks
  data = read.csv(infile) 
  info = describe(data)
  cnames = colnames(data)
    
  ## if user has specified the columns to be used
  if (!is.null(cols)){
    verifyCols <- (cols %in% cnames)
    if (!all(verifyCols)) stop(paste0('>', cols[!verifyCols], '< is not in the data columns.', collapse = '\n'))
    tmpData <- data[,cols]
  
  ## try to geuss the columns that is be used
  } else {
    info$cols.bait <- grepl(paste(bait, collapse='.*'), cnames) & (!info$cols.ratios)
    dataBait <- data[,info$cols.bait]
    dataMock <- data[,info$cols.control]
    dataComb <- cbind(dataBait,dataMock)
    dataComb <- dataComb[,c(1,3,2,4)] # iTRAQ duplicates only for now // numbers should be checked
    tmpData <- cbind(data[info$col.accession], dataComb)
    
  }
  
  # Replace zeros with NAs
  tmpData[tmpData == 0] <- NA
  info$count.na <- sum(as.numeric(is.na(tmpData)))
  
  # 1) transform the data and 2) median normalizeation
  tmpData = pTransform(tmpData, type = transform)
  tmpData = normalize(tmpData, type = normalization)
  
  # 3) remove non human proteins and proteins with < 2 unique peptides
  rowsEnoughProteins <- !(data[,info$col.unique.proteins] >= 2) 
  tmpData <- tmpData[!rowsEnoughProteins,]
  info$rows.with.too.few.proteins <- sum(as.numeric(rowsEnoughProteins))
  rowsNotHuman <- !grepl(filter, tmpData[,1])
  tmpData <- tmpData[!rowsNotHuman,]
  info$rows.removed.by.filter <- sum(as.numeric(rowsNotHuman))
      
  # 4) extract gene only
  tmpData[,1] <- strSplitGene(tmpData[,1])
      
  # 5) impute if needed
  tmpData = tmpData[complete.cases(tmpData),]
  info$total.rows.remove <- nrow(data) - nrow(tmpData)
  
  # 6) calculate log fold change
  tmpData = logFC(tmpData)
  
  ## resulting data
  if (!verbose){
    return(tmpData)
  } else {
    return(list(data=tmpData, info=info))
  }

}




