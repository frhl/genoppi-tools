#' @title pipeline for process proteomics data
#' @description a function that will process ip data
#' from different baits, files and cells
#' @param bait a vector containing that contains the bait and cell type that should be matched in a column.
#' @param infile the file path or data.frame that contains the raw data, i.e. accession numbers, intensity values, ratios etc.
#' @param cols optional manual entry. A vector of columns that are present in the dataset header. Follows 
#' the format of columns: acession, bait1, mock1, bait2, mock2, bait3, ..
#' @param impute how should missing data be imputed? NULL means that missing data rows are dropped. Will take a list
#' containing shift and stdwidth for gaussian imputation. For further details, see \code{?impute}. 
#' @param transform charcacter. an R-command for how the data should be transformed.
#' @param normalization character. an R-command for how the data should be collumn-wise transformned.
#' @param filter character. only accession IDs of the filter specified are included.
#' @param raw will return the data.frame alongside the raw intensity values.
#' @param verbose boolean. if true, returns the table and a list with updates.
#' @return a table that can be inputted to genoppi

prepare <- function(bait, infile, cols = NULL, impute = list(stdwidth = 0.5, shift = -0.8), 
                    transform = 'log2', normalization = 'median', filter = "HUMAN", raw = F, verbose = F){
  
  ## do some initial checks
  if (is.character(infile)) data = read.csv(infile) else data = as.data.frame(infile)
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
    stopifnot(!is.null(dataBait) | !is.null(dataMock)) # should not be null
    stopifnot(ncol(dataBait) == ncol(dataMock)) # should have same amount of columns
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
  
  ## note (this is ugly code.. clean up)
  tmpData$enoughProteins <- data[,info$col.unique.proteins] >= 2
  tmpData <- tmpData[tmpData$enoughProteins == TRUE,]
  info$rows.with.too.few.proteins <- sum(as.numeric(!tmpData$enoughProteins))
  tmpData$human <- grepl(filter, tmpData[,1])
  tmpData <- tmpData[tmpData$human,]
  info$rows.removed.by.filter <- sum(as.numeric(!tmpData$human))
      
  # 4) extract gene only
  tmpData[,1] <- strSplitGene(tmpData[,1])
      
  # 5) impute if needed
  if (is.null(impute)) {
    tmpData = tmpData[complete.cases(tmpData),] 
  } else {
      if (all(c('stdwidth' , 'shift') %in% names(impute))){
        tmpData = impute.gaussian(tmpData, impute$stdwidth, impute$shift)
      } else {stop('Use impute params "stdwidth" and "shift" only.')}
    }

  # 6) calculate log fold change
  info$total.rows.remove <- nrow(data) - nrow(tmpData)
  tmpData = logFC(tmpData)
  
  # clean out the data and remove intensity columns
  if (!raw){
    tmpData = tmpData[,grepl('rep|Acc|impute', colnames(tmpData))]
  }
  
  ## resulting data
  if (!verbose){
    return(tmpData)
  } else {
    return(list(data=tmpData, info=info))
  }

}




