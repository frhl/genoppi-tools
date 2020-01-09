## initial analysis
if (F){
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
(files = list.files('~/Projects/03_MICOM/data/raw/martin10/', full.names = T))

########################################################### 
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: FN1 (recovered as alias:FINC)
# cell: SMC
# note: three different isoforms are detected, but they all have the exact same intensity
# values.
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin10//10martin_SMC_FN1_proteins.csv"
head(read.csv(infile))
bait = 'FN' #FN1
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'CIG|MSF|FNZ|LETS|GFND|FINC') # Alias is FINC...
read.csv(infile)[det,]  # three isoforms.. But they all have some intensity values.
dirs <- mkdir(bait = bait, cell = cell, run = 'martin10') 

# data processing
data <- prepare(c(cell, bait), infile = infile, verbose = T, filter.ignore = 'FINC') %>% mttest()
data$gene <- as.character(data$gene)
data[data$gene == '!FINC', ]$gene <- 'FN1'

write.table(data,file=dirs$txtpath,row.names=F,sep="\t",quote=F)
pdf(dirs$pdfpath, height=8, width=8)

## plot volcano and scatter plot
data %>% designate(FDR < 0.1) %>% plotVolcano('FN1', title = 'FN1 vs Control (EC).')
r = data %>% designate(FDR < 0.1) %>% plotScatter('FN1', title = 'FN1 (EC) Replicate Correlation')

### plot inweb interactors and proxies
known.interactors = interactors('FN1', T)
data %>% designate(FDR < 0.1) %>% plotOverlap('FN1', known.interactors, 'FN1 vs Control (EC): InWeb Overlap')
data %>% designate(FDR < 0.1, logFC > 0) %>% plotOverlap('FN1', known.interactors, 'FN1 vs Control (EC): InWeb Overlap')

# plot proxies
proxies <- expandProxies(data)
data %>% designate(FDR < 0.1, logFC > 0) %>% plotProxies(proxies)
graphics.off()

# save overview
report = list(bait='FN1', 
              cell=cell, 
              corr = r, 
              file=infile, 
              data = data, 
              data.bait.found = TRUE,
              data.bait.enrichment = data[data$gene == 'FN1', ]$logFC[1],
              interactors = known.interactors, 
              in.data.interactors = get.indata.interactors(data, known.interactors),
              significant.interactors = get.interactors.genes(data, known.interactors))

write.csv(review(report), file = dirs$csvpathsum)  


########################################################### 
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: JCAD (recovered)
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin10//10martin_SMC_JCAD_proteins.csv" 
head(read.csv(infile))
bait = 'JCAD'
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'JCAD|KIAA1462') 
read.csv(infile)[det,] 
dirs <- mkdir(bait = bait, cell = cell, run = 'martin10') 

# data processing
data <- prepare(c(cell, bait), infile = infile) %>% mttest()
write.table(data,file=dirs$txtpath,row.names=F,sep="\t",quote=F)
pdf(dirs$pdfpath, height=8, width=8)

## plot volcano and scatter plot
data %>% designate(FDR < 0.1) %>% plotVolcano('JCAD', title = 'JCAD vs Control (SMC).')
r = data %>% designate(FDR < 0.1) %>% plotScatter('JCAD', title = 'JCAD (SMC) Replicate Correlation')

### plot inweb interactors and proxies
known.interactors = interactors('KIAA1462', T) # alias for JCAD: https://www.genecards.org/cgi-bin/carddisp.pl?gene=JCAD
data %>% designate(FDR < 0.1) %>% plotOverlap('JCAD', known.interactors, 'JCAD vs Control (EC): InWeb Overlap')
data %>% designate(FDR < 0.1, logFC > 0) %>% plotOverlap('JCAD', known.interactors, 'JCAD vs Control (EC): InWeb Overlap')

# plot proxies
proxies <- expandProxies(data)
data %>% designate(FDR < 0.1, logFC > 0) %>% plotProxies(proxies)
graphics.off()

# save overview
report = list(bait='JCAD (alias:KIAA1462)', 
              cell=cell, 
              corr = r, 
              file=infile, 
              data = data, 
              data.bait.found = FALSE,
              data.bait.enrichment = NA,
              interactors = known.interactors, 
              in.data.interactors = get.indata.interactors(data, known.interactors),
              significant.interactors = get.interactors.genes(data, known.interactors))

write.csv(review(report), file = dirs$csvpathsum)

                    
                    