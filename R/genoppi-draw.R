rm(list=ls())
library(limma)
library(ggplot2)
library(ggrepel)
library(hash)

# load Genoppi functions
source("/psych/lagelab/yuhanhsu/Genoppi/src/genoppiFunctions.r")

# input arguments
args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1] # Genoppi input table
fc_cutoff <- args[2] # "positive", "negative", or "both"
p_cutoff <- args[3] # numeric value or "NA"
fdr_cutoff <- args[4] # numeric value or "NA"
bait_name <- args[5] # bait name
imp_list_file <- args[6]
gene_lists_file <- args[7] # file containing gene list file names, or "NA"
output_stats_file <- args[8] # output file containing moderated t-test stats
output_plots_file <- args[9] # output file containing plots

#input_file <- 'GenoppiInput/BCL2.HDFNvsG401.NoImp.GenoppiInput.txt' 
#fc_cutoff <- 'both'
#p_cutoff <- 'NA'
#fdr_cutoff <- 0.1
#bait_name <- 'BCL2'
#imp_list_file <- 'NA'
#gene_lists_file <- 'test.GenoppiIP.genelists.txt'
#output_stats_file <- 'test.BCL2_HDFNvsG401.GenoppiOutput.txt'
#output_plots_file <- 'test.BCL2_HDFNvsG401.GenoppiPlots.pdf'

# read in GenoppiInput file: gene name followed by logFC for each replicate
input_df <- read.table(input_file,header=T,sep="\t")

# calculate and output moderated t-test statistics
out_stats <- calculate_moderated_ttest(input_df)
write.table(out_stats,file=output_stats_file,row.names=F,sep="\t",quote=F)

# define signifcant proteins
out_stats$significant <- TRUE # fc_cutoff=="both"
if (fc_cutoff=="positive") { 
  out_stats$significant <- out_stats$logFC > 0
} else if (fc_cutoff=="negative") { 
  out_stats$significant <- out_stats$logFC < 0 }

significant_pos <- out_stats$logFC > 0
significant_neg <- out_stats$logFC < 0

if (p_cutoff!="NA") {
  p_cutoff <- as.numeric(p_cutoff)
  out_stats$significant <- out_stats$significant & out_stats$pvalue < p_cutoff
  significant_pos <- significant_pos & out_stats$pvalue < p_cutoff
  significant_neg <- significant_neg & out_stats$pvalue < p_cutoff
}
if (fdr_cutoff!="NA") {
  fdr_cutoff <- as.numeric(fdr_cutoff)
  out_stats$significant <- out_stats$significant & out_stats$FDR < fdr_cutoff
  significant_pos <- significant_pos & out_stats$FDR < fdr_cutoff
  significant_neg <- significant_neg & out_stats$FDR < fdr_cutoff	
}

# generate plots and run enrichment tests
pdf(output_plots_file,height=4,width=4)

draw_volcano(out_stats, bait_name)
draw_scatter(out_stats, bait_name)


# label IMPUTED genes on volcano plot
if (imp_list_file != 'NA') {
  impDf <- read.table(imp_list_file,header=T,sep="\t")
  draw_overlap(out_stats,bait_name,'IMPUTED',impDf,FALSE)
}

# InWeb overlap enrichment (only run if bait in InWeb)
load("../BioID_benchmark/InWeb_combined_Oct2018.RData") # inweb_hash hash object

if (bait_name %in% keys(inweb_hash)) {
  inwebDf <- data.frame(gene=keys(inweb_hash))
  inwebDf$significant <- inwebDf$gene %in% inweb_hash[[bait_name]]
  
  # significance direction "both"
  draw_overlap(out_stats,bait_name,'InWeb',inwebDf,TRUE)
  draw_overlap(out_stats,bait_name,'InWeb',inwebDf,FALSE)
  
  # significance direction "positive"
  temp_stats <- out_stats
  temp_stats$significant <- significant_pos
  draw_overlap(temp_stats,bait_name,'InWeb, +',inwebDf,TRUE)
  draw_overlap(temp_stats,bait_name,'InWeb, +',inwebDf,FALSE)
  
  # significance direction "negative"
  temp_stats$significant <- significant_neg
  draw_overlap(temp_stats,bait_name,'InWeb, -',inwebDf,TRUE)
  draw_overlap(temp_stats,bait_name,'InWeb, -',inwebDf,FALSE)
  
}

# with gene lists overlap enrichment
if (gene_lists_file != 'NA') {
  geneListTable <- read.table(gene_lists_file,header=T,sep="\t",stringsAsFactors=F)
  # loop through each gene list
  for (i in 1:dim(geneListTable)[1]) {
    listDf <- read.table(geneListTable[i,2],header=T,sep="\t")
    draw_overlap(out_stats,bait_name,geneListTable[i,1],listDf,TRUE)
    
    # significance direction "positive"
    temp_stats$significant <- significant_pos
    draw_overlap(temp_stats,bait_name,paste(geneListTable[i,1],', +',sep=''),listDf,TRUE)	
    
    # significance direction "negative"
    temp_stats$significant <- significant_neg
    draw_overlap(temp_stats,bait_name,paste(geneListTable[i,1],', -',sep=''),listDf,TRUE)
    
  }
}

dev.off()
