#' @title Update data.bases
#' @description A function for generating a hash table that contains the columns 'primary' and 'synonym'
#' 
#' 


mkhashtable <- function(dat, split = ' ', chars.remove = '\\*'){
  
  require(hashmap)
  write('making hashmap..',stdout())
  
  # check input and order it
  stopifnot('primary' %in% colnames(dat))
  stopifnot('synonym' %in% colnames(dat))
  dat = dat[,c('primary', 'synonym')]
  tmp <- c(NA,NA)
  
  # iterate through all data sequentially
  for (i in 1:nrow(dat)){
    raw.primary = as.character(gsub(chars.remove, '', dat$primary[i]))
    raw.synonym = as.character(dat$synonym[i])

    # we assume that there is always at least 
    # one letter in the synonym or primary ID
    if (grepl('[A-Za-z]', raw.primary)){
      if (grepl('[A-Za-z]', raw.synonym)){
        split.synonym = unique(unlist(strsplit(raw.synonym, split)))
        for (j in split.synonym){
          splitted <- gsub(chars.remove,'',j)
          splitted <- gsub(' ','',splitted)
          tmp <- rbind(tmp, c(splitted, raw.primary))
        }
      } else {
        tmp <- rbind(tmp, c(raw.primary, raw.primary))
      }
    }
  }
  
  # modify data.frame and hash it
  tmp <- as.data.frame(tmp[complete.cases(tmp),])
  colnames(tmp) <- c('alias', 'primary')
  rownames(tmp) <- NULL
  #return(tmp)
  return(hashmap(keys=as.character(tmp$alias), values=as.character(tmp$primary)))
}


## make uniprot hashtable
#dat <- as.data.frame(fread('~/Desktop/geneMapping/uniprot-reviewed_yes+AND+proteome_up000005640.tab'))
#dat <- dat[,colnames(dat) %in% c('Gene names  (primary )', 'Gene names  (synonym )')]
#dat <- dat[,c(2,1)]
#colnames(dat) <- c('primary', 'synonym')
#uniprot_aliases <- mkhashtable(dat)
#save(uniprot_aliases, file = '~/Toolbox/packages/pRoteomics/data/uniprot_aliases.RData')
#hm <- hashmap(keys=as.character(db$alias), values=as.character(db$protein))

#.libPaths('~/Toolbox/rlib/')
#library('biomaRt')
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#filters = listFilters(ensembl)
#dat = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'uniprot_gn_symbol'), mart = ensembl)
#dat = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'hgnc_id', 'uniprot_gn_symbol', 'uniprot_gn_id'), mart = ensembl)
#dat1 <- dat

#dat2 <- data.frame(primary=dat1$ensembl_gene_id, synonym=paste(dat1$hgnc_symbol, dat1$uniprot_gn_symbol, sep = ' ') )
#to_ensembl <- mkhashtable(dat2)
#hashmap(keys=as.character(tmp$alias), values=as.character(tmp$protein))


