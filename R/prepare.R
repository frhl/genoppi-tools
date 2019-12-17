#' @title pipeline for preparing genomics data for statistical analysis
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
#' @param firstcol will change the name of the first column to the string indicated
#' @param filter.ignore will try to match the inputted vector or character to acession IDs. If sucessful,
#' it will ignore further filtering of this item. This could for instance be used, if the bait only has
#' one unique protein, and would therefore otherwise be filtered.
#' @param verbose boolean. if true, returns the table and a list with updates.
#' @export
#' @return a table that can be inputted to genoppi

prepare <- function(bait, infile, cols = NULL, impute = list(stdwidth = 0.5, shift = -1.8), 
                    transform = 'log2', normalization = 'median', filter = "HUMAN", raw = F, firstcol = 'gene',
                    filter.ignore = NULL, verbose = F){
  
  ## do some initial checks
  ## and read in the file
  if (all(is.null(bait))) stop('Bait can not be NULL!')
  #if (is.null(dim(infile)) & !file.exists(infile)) stop('Infile must be either a file path that exists or a data.frame.')
  if (is.character(infile)) data = read.csv(infile) else data = as.data.frame(infile)
  info = describe(data)
  cnames = colnames(data)
      
  ## if user has specified the columns to be used
  if (!is.null(cols)){

    verifyCols <- (cols %in% cnames)
    if (!all(verifyCols)) stop(paste0('>', cols[!verifyCols], '< is not in the data columns.', collapse = '\n'))
    tmpData <- data[,cols]
  
  } else {
    ## try to geuss the columns that is be used
    baitFound <- !unlist(lapply(bait, function(x) any(grepl(x, cnames))))
    info$cols.bait <- grepl(paste(bait, collapse='.*'), cnames) & (!info$cols.ratios) & (!info$cols.control)
    dataBait <- data[,info$cols.bait]
    dataMock <- data[,info$cols.control]
    
    # check format and give transparent error message
    if (any(baitFound)) stop(paste(c(bait[baitFound], '(bait) not in data columns!.'), collapse = ' '))
    if (sum(info$cols.bait) == 1) stop('expected at least two columns of baits, only one was found!')
    if (sum(info$cols.mock) == 1) stop('expected at least two columns of controls, only one was found!')
    if (!ncol(dataBait)) stop('bait columns were not found!')
    if (!ncol(dataMock)) stop('mock columns were not found!')
    if (is.null(dataBait) | is.null(dataMock)) stop('disproprionate amount of bait and mock columns were found')
    if (ncol(dataBait) != ncol(dataMock)) stop('disproprionate amount of bait and mock columns were found')
    
    # verbose
    if (verbose) warn(paste('[Verbose] Selected bait cols:', paste(cnames[info$cols.bait], collapse = ' ')))
    if (verbose) warn(paste('[Verbose] Selected mock cols:', paste(cnames[info$cols.control], collapse = ' ')))
    
    # prepare data # should have same amount of columns
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
  tmpData$filter.ignore <- detect(tmpData, filter.ignore) # allow user to ignore rows
  nignored <- sum(tmpData$filter.ignore)
  tmpData$enoughProteins <- data[,info$col.unique.proteins] >= 2 | tmpData$filter.ignore
  tmpData <- tmpData[tmpData$enoughProteins == TRUE,]
  tmpData$human <- grepl(filter, tmpData[,1]) | tmpData$filter.ignore
  tmpData <- tmpData[tmpData$human,]
  if (verbose & nignored > 0) warn(paste('[filtering]', nignored, 'entries was ignored.'))
  
  # 4) convert from uniprot to HGNC\
  matr <- acession.matrix(tmpData[,1]) # first column is the acession
  matr.convert <- acession.convert(matr, verbose = verbose)
  tmpData$Accession <-matr.convert$hgnc # extract hgnc symbol
  
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
    if (!is.null(firstcol)) colnames(tmpData)[1] <- 'gene'
    return(tmpData)
  } else (
    return(list(data=tmpData, info=info))
  )
  
}




