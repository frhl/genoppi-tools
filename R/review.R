

review <- function(params = NULL){
  
  lst <- list()
  
  if (!is.null(params)){
    lstnames <- names(params)
    
    if ('experiment' %in% lstnames){
      lst$experiment <- params$experiment
    }
    if ('bait' %in% lstnames){
      lst$bait <- params$bait
    }
    if ('cell' %in% lstnames){
      lst$bait <- params$cell
    }
    # How many of uniprots were converted sucessfully to HGNC
    if ('coverage' %in% lstnames){
      lst$bait <- params$coverage
    }
    
    ## sumamrize data itself
    if ('data' %in% lstnames){
      lst$data.rows <- nrow(params$data) # proteins indentified
      lst$data.imputed <- sum(as.numeric(params$data$imputed)) # data imputed
      lst$data.fdr.ok <- sum(as.numeric(params$data$FDR < 0.1)) # data with FDR < 0.1
      lst$data.nominal.sig <- sum(as.numeric(params$data$FDR < 0.05)) # data with pvalue < 0.05
      lst$data.pos.logfc <- sum(as.numeric(params$data$logFC > 0)) # data with positive logFC
      lst$data.pos.log.fc.and.fdr.ok <- sum(as.numeric(params$data$logFC > 0 & params$data$FDR < 0.1)) # data with both positive logFC and FDR < 0.1
      lst$data.bait.found <- params$bait %in% params$data$gene # can the bait be recovered?
      if (lst$data.bait.found) lst$data.bait.enrichment <- params$data[params$data$gene == bait,]$logFC else lst$data.bait.found <- NA # bait enrichment (inlogFC)
      #data.rep.corr <- cor(params$data$rep1, params$data$rep2)
    }
    ## summarize inweb interactors
    if ('interactors' %in% lstnames & 'data' %in% lstnames){
      lst$interactors.n <- nrow(params$interactors[params$interactors$significant == TRUE, ]) # inweb interactors found
      lst$interactors.in.data <- sum(as.numeric(params$interactors[params$interactors$significant == TRUE, ]$gene %in% params$data$gene)) # inweb interactors actually in data
    }
    if ('corr' %in% lstnames){
      lst$rep.cor <- params$corr$r
    }
    if ('infile' %in% lstnames){
      lst$file <- params$infile
    }

    
  }
  
  lst <- do.call(data.frame, lst)
  return(lst)
  
  
  
  
  
}