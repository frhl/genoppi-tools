
if (F){
  ## initial analysis
  .libPaths('~/Toolbox/rlib')
  setwd('~/Projects/03_MICOM/')
  library(dplyr)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(hash)
  library(dplyr)
  
  #library(rProteomics, lib.loc = '~/Toolbox/rlib/')
  devtools::load_all('~/Toolbox/packages/pRoteomics/')
}
(files = list.files('~/Projects/03_MICOM/data/raw/martin6/', full.names = T))

########################################################### 
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: HDAC9
# control: control
# cell: SMC
###########################################################


infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin6//6martin_HDAC9_proteins.csv" 
head(read.csv(infile))
bait = 'HDAC9'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin6') 

# chcek if bait in data
det <- detect(infile, bait); any(det)
read.csv(infile)[det, ]

# Prepare data for analysis by getting replicate folds
data <- prepare(c(cell, bait), infile = infile) %>% mttest()
write.table(data,file=dirs$txtpath,row.names=F,sep="\t",quote=F)
pdf(dirs$pdfpath, height=8, width=8)

## plot volcano and scatter plot
data %>% designate(FDR < 0.1) %>% plotVolcano(bait, title = 'HDAC9 vs Control (SMC). Note: Only one unique peptide detected.')
r = data %>% designate(FDR < 0.1) %>% plotScatter(bait, title = 'HDAC9 Replicate Correlation')

# plot imputed data
#data %>% designate(imputed == TRUE, FDR < 0.1) %>% plotVolcano(bait, title = 'Imputation summary', sub2 ='imputed')

### plot inweb interactors
known.interactors = interactors('HDAC9', T)
data %>% designate(FDR < 0.1) %>% plotOverlap(bait, known.interactors, 'HDAC9 vs Control (EC): InWeb Overlap')

## find proxy baits
proxies <- expandProxies(data)
data %>% designate(FDR < 0.1) %>% plotProxies(proxies)
graphics.off()

# save overview
report = list(bait=bait, 
              cell=cell, 
              corr = r, 
              file=infile, 
              data = data, 
              data.bait.found = TRUE,
              data.bait.enrichment = data[data$gene == bait, ]$logFC,
              data.bait.enriched = data[data$gene == bait, ]$FDR < 0.1 & data[data$gene == bait, ]$logFC > 0,
              interactors = known.interactors, 
              in.data.interactors = get.indata.interactors(data, known.interactors),
              significant.interactors = get.interactors.genes(data, known.interactors))

write.csv(review(report), file = dirs$csvpathsum)