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
#' @family genoppi
#' @export


genoppi <- function(input_file,  bait_name, output_stats_file, output_plots_file, gene_lists_file = NA,
                    fc_cutoff = 'both', fdr_cutoff = 0.1, p_cutoff = NA, imp_list_file = NA, verbose = T, 
                    aggregate_info = T, debug = F){

  require(limma)
  require(ggplot2)
  require(ggrepel)
  require(hash)
  
  if (debug) browser()
  
  # initial checks
  stopifnot(length(bait_name) == 1)
  
  # read in GenoppiInput file: gene name followed by logFC for each replicate

  input_df <- read.table(input_file,header=T,sep="\t")
  stopifnot(any(grepl('rep', colnames(input_df))))
  
  # calculate and output moderated t-test statistics
  out_stats <- genoppi.ttest(input_df, keep = c('imputed'))
  write.table(out_stats,file=output_stats_file,row.names=F,sep="\t",quote=F)
  
  # define significant proteins
  out_stats$significant <- TRUE # fc_cutoff=="both"
  if (fc_cutoff=="positive") { 
  	out_stats$significant <- out_stats$logFC > 0
  } else if (fc_cutoff=="negative") { 
  	out_stats$significant <- out_stats$logFC < 0 }
  
  significant_pos <- out_stats$logFC > 0
  significant_neg <- out_stats$logFC < 0
  
  #browser()
  
  ## check whether any preys are significant by p_value
  if (!is.na(p_cutoff)) {
  	p_cutoff <- as.numeric(p_cutoff)
  	out_stats$significant <- out_stats$significant & out_stats$pvalue < p_cutoff
  	significant_pos <- significant_pos & out_stats$pvalue < p_cutoff
  	significant_neg <- significant_neg & out_stats$pvalue < p_cutoff
  	if (verbose){
  	  p_cutoff_n_targets = length(out_stats$significant)
  	  p_cutoff_n_sig = sum(as.numeric(out_stats$significant))
  	  write(paste('[PVAL]:', p_cutoff_n_sig, 'of', p_cutoff_n_targets,'target(s) are significant with P-value <', p_cutoff),stdout()) 
  	}
  }

  ## check whether any preys are significant by FDR
  if (!is.na(fdr_cutoff)) {
  	fdr_cutoff <- as.numeric(fdr_cutoff)
  	out_stats$significant <- out_stats$significant & out_stats$FDR < fdr_cutoff
  	significant_pos <- significant_pos & out_stats$FDR < fdr_cutoff
  	significant_neg <- significant_neg & out_stats$FDR < fdr_cutoff
  	if (verbose) {
  	  fdr_cutoff_n_targets = length(out_stats$significant)
  	  fdr_cutoff_n_sig = sum(as.numeric(out_stats$significant))
  	  write(paste('[FDR]:', fdr_cutoff_n_sig, 'of', fdr_cutoff_n_targets,'target(s) are significant with FDR <', fdr_cutoff),stdout())
  	}
  }
  
  # generate plots and run enrichment tests
  pdf(output_plots_file,height=4,width=4)
  plotVolcano(out_stats, bait_name)
  plotScatter(out_stats, bait_name)
  
  # label IMPUTED genes on volcano plot
  #if (!is.na(imp_list_file)) {
  	#impDf <- read.table(imp_list_file,header=T,sep="\t")
  
  
  ## this function must be dissected
  #browser()
  impDf <- out_stats[out_stats$impute == TRUE, ]
  plotOverlap(out_stats, bait_name, 'imputed' ,impDf ,TRUE, debug = F)
  if (verbose) write(paste('[IMPUTE]:', nrow(impDf), 'of', nrow(out_stats),'target(s) are significant and imputed.'),stdout()) 
  
  
  # InWeb overlap enrichment (only run if bait in InWeb)
  data(inweb_hash)
  
  if (bait_name %in% keys(inweb_hash)) {
    #browser()
    
  	inwebDf <- data.frame(gene=keys(inweb_hash))
  	inwebDf$significant <- inwebDf$gene %in% inweb_hash[[bait_name]]
  
  	# significance direction "both"
  	plotOverlap(out_stats,bait_name,'InWeb',inwebDf,TRUE, debug = F)
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



  



