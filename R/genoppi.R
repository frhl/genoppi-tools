#' @title Genoppi analysis pipeline
#' @description the genoppi combined analysis script. This will run the pipeline
#' generating correlation plots, 
#' @param input_file Genoppi input table
#' @param fc_cutoff character "positive", "negative", or "both"
#' @param p_cutoff numeric value or "NA"
#' @param fdr_cutoff numeric value or "NA"
#' @param bait_name bait name
#' @param gene_lists_file file containing gene list file names, or "NA"
#' @param imp_list_file file containing an overview over which proteins were imputed
#' @param output_stats_file output file containing moderated t-test stats (must be a .txt file)
#' @param output_plots_file output file containing plots (must be a .pdf file)
#' @author aprilkim/Yuhanshu/flassen
#' @note functions modified from Genoppi source code
#' @export

genoppi <- function(input_file,  bait_name, output_stats_file, output_plots_file, gene_lists_file = NA,
                    fc_cutoff = 'both', fdr_cutoff = 0.1, p_cutoff = NA, imp_list_file = NA, debug = F){

  require(limma)
  require(ggplot2)
  require(ggrepel)
  require(hash)
  
  if (debug) browser()
  
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
  stopifnot(any(grepl('rep', colnames(input_df))))
  
  # calculate and output moderated t-test statistics
  out_stats <- moderatedTTest(input_df)
  write.table(out_stats,file=output_stats_file,row.names=F,sep="\t",quote=F)
  
  # define signifcant proteins
  out_stats$significant <- TRUE # fc_cutoff=="both"
  if (fc_cutoff=="positive") { 
  	out_stats$significant <- out_stats$logFC > 0
  } else if (fc_cutoff=="negative") { 
  	out_stats$significant <- out_stats$logFC < 0 }
  
  significant_pos <- out_stats$logFC > 0
  significant_neg <- out_stats$logFC < 0
  
  if (!is.na(p_cutoff)) {
  	p_cutoff <- as.numeric(p_cutoff)
  	out_stats$significant <- out_stats$significant & out_stats$pvalue < p_cutoff
  	significant_pos <- significant_pos & out_stats$pvalue < p_cutoff
  	significant_neg <- significant_neg & out_stats$pvalue < p_cutoff
  }
  if (!is.na(fdr_cutoff)) {
  	fdr_cutoff <- as.numeric(fdr_cutoff)
  	out_stats$significant <- out_stats$significant & out_stats$FDR < fdr_cutoff
  	significant_pos <- significant_pos & out_stats$FDR < fdr_cutoff
  	significant_neg <- significant_neg & out_stats$FDR < fdr_cutoff	
  }
  
  # generate plots and run enrichment tests
  pdf(output_plots_file,height=4,width=4)
  plotVolcano(out_stats, bait_name)
  plotScatter(out_stats, bait_name)
  
  # label IMPUTED genes on volcano plot
  #if (!is.na(imp_list_file)) {
  	#impDf <- read.table(imp_list_file,header=T,sep="\t")
  
  
  ## this function must be dissected
  impDf <- out_statsß
  plotOverlap(out_stats,bait_name,'imputed',impDf,FALSE, debug = T)
  #}
  
  # InWeb overlap enrichment (only run if bait in InWeb)
  data(inweb_hash)
  
  if (bait_name %in% keys(inweb_hash)) {
    
  	inwebDf <- data.frame(gene=keys(inweb_hash))
  	inwebDf$significant <- inwebDf$gene %in% inweb_hash[[bait_name]]
  
  	# significance direction "both"
  	plotOverlap(out_stats,bait_name,'InWeb',inwebDf,TRUE)
  	plotOverlap(out_stats,bait_name,'InWeb',inwebDf,FALSE)
  
  	# significance direction "positive"
  	temp_stats <- out_stats
  	temp_stats$significant <- significant_pos
  	plotOverlap(temp_stats,bait_name,'InWeb, +',inwebDf,TRUE)
  	plotOverlap(temp_stats,bait_name,'InWeb, +',inwebDf,FALSE)
  	
  	# significance direction "negative"
  	temp_stats$significant <- significant_neg
  	plotOverlap(temp_stats,bait_name,'InWeb, -',inwebDf,TRUE)
  	plotOverlap(temp_stats,bait_name,'InWeb, -',inwebDf,FALSE)
  	
  }
  
  # with gene lists overlap enrichment
  if (!is.na(gene_lists_file)) {
  	geneListTable <- read.table(gene_lists_file,header=T,sep="\t",stringsAsFactors=F)
  	
  	# loop through each gene list
  	for (i in 1:dim(geneListTable)[1]) {
  		listDf <- read.table(geneListTable[i,2],header=T,sep="\t")
  		plotOverlap(out_stats,bait_name,geneListTable[i,1],listDf,TRUE)
  	
  		# significance direction "positive"
  		temp_stats$significant <- significant_pos
  		plotOverlap(temp_stats,bait_name,paste(geneListTable[i,1],', +',sep=''),listDf,TRUE)	
  	
  		# significance direction "negative"
  		temp_stats$significant <- significant_neg
  		plotOverlap(temp_stats,bait_name,paste(geneListTable[i,1],', -',sep=''),listDf,TRUE)
  	
  	}
  }
  
  graphics.off()
}
