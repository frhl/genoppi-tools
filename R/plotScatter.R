#' @title Draw Replicate Scatter Plot
#' @description Draws replicate correlation scatter plot(s)
#' @param df something
#' @param bait something
#' @author April Kim / Frederik Heymann
#' @family genoppi
#' @export


### --------------------------------------------------
### replicate correlation scatter plot(s)
plotScatter <- function(df, bait, title = ''){
  
  # iterate through each pair of replicates
  #browser()

  #nRep <- dim(df)[2] - 5 # minus gene, logFC, pvalue, FDR, significant columns
  colRep <- as.vector(grepl('rep', colnames(df)) & unlist(lapply(df, is.numeric)))
  nRep <- sum(as.numeric(colRep))
  
  for (i in 1:(nRep-1)) {
    for (j in (i+1):nRep) {
      
      col1 <- paste("rep",i,sep="")
      col2 <- paste("rep",j,sep="")
      r <- cor(df[,col1],df[,col2]) # Pearson correlation
      temp_df <- data.frame(gene=df$gene,temp_rep1=df[,col1],temp_rep2=df[,col2],significant=df$significant)
      
      # start scatterplot
      p <- ggplot(temp_df, aes(x=temp_rep1, y=temp_rep2)) +
        
        # plot all proteins (green = significant, blue = not significant)
        geom_point(alpha=0.5, size=1.5, color=ifelse(df$significant, "springgreen3", "royalblue2")) +
        
        # label bait (red = signficant, orange = not significant)
        geom_point(subset(temp_df, gene==bait & significant), mapping=aes(x=temp_rep1, y=temp_rep2),
                   size=2, color="red") + 
        geom_point(subset(temp_df, gene==bait & !significant), mapping=aes(x=temp_rep1, y=temp_rep2),
                   size=2, color="orange") +
        geom_point(subset(temp_df, gene==bait), mapping=aes(x=temp_rep1, y=temp_rep2),
                   size=2, color="black", shape=1) +	
        geom_text_repel(subset(temp_df, gene==bait), mapping=aes(label=gene),
                        arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.15, "lines"),
                        point.padding=unit(0.2, "lines"), color="black", size=3) +
        
        # identity line, title (with correlation), theme
        geom_abline(intercept=0, slope=1, linetype="longdash", size=0.2) +
        labs(title = title, subtitle = paste("correlation:",format(r,digits=3))) + xlab(col1) + ylab(col2) +
        theme_bw() + theme(axis.line=element_line(color="grey")) + ggstamp()
      
      print(p)
      
    }
  }
  return(list(r=r))
}
