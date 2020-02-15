#' @title Draw a scatter plot of replicates
#' @description Draws replicate correlation scatter plot(s)
#' @param df A data.frame containing replicate columns.
#' @param bait The bait that was used. A string.
#' @param title Add a title to the ggplot.
#' @author April Kim / Frederik Heymann
#' @family genoppi
#' @export

plotScatterMinimal <- function(df, bait, title = '', size_point = 3, size_text=3, color_alpha=0.8){
  
  # inital checks of input
  stopifnot(any(grepl('rep', colnames(df))))
  colRep <- as.vector(grepl('rep', colnames(df)) & unlist(lapply(df, is.numeric)))
  reps = regmatches(colnames(df), regexpr('rep[0-9]',colnames(df)))
  nRep <- sum(as.numeric(colRep))
  
  require(ggplot2)
  require(ggrepel)
  
  for (i in 1:(nRep-1)) {
    for (j in (i+1):nRep) {
      
      # set columns of replicates
      col1 <- reps[i]#paste("rep",i,sep="")
      col2 <- reps[j]#paste("rep",j,sep="")
      r <- cor(df[,col1],df[,col2]) # Pearson correlation
      temp_df <- data.frame(gene=df$gene,temp_rep1=df[,col1],temp_rep2=df[,col2],significant=df$significant)
      
      # start scatterplot
      p <- ggplot(temp_df, aes(x=temp_rep1, y=temp_rep2)) +
        
        # plot all proteins (green = significant, blue = not significant)
        geom_point(alpha=0.95, size=size_point+0.3, color=ifelse(df$significant, "black", "grey"), shape=ifelse(df$significant, 1, 1)) +
        geom_point(alpha=color_alpha, size=size_point, color=ifelse(df$significant, "springgreen3", "grey")) +
        
        
        # label bait (red = signficant, orange = not significant)
        geom_point(subset(temp_df, gene==bait & significant), mapping=aes(x=temp_rep1, y=temp_rep2),size=size_point, color="red") + 
        geom_point(subset(temp_df, gene==bait & !significant), mapping=aes(x=temp_rep1, y=temp_rep2),size=size_point, color="orange") +
        geom_point(subset(temp_df, gene==bait), mapping=aes(x=temp_rep1, y=temp_rep2), size=size_point, color="black", shape=1) +	
        geom_text_repel(subset(temp_df, gene==bait), mapping=aes(label=gene),
                        arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.15, "lines"),
                        point.padding=unit(0.2, "lines"), color="black", size=size_text) +
        
        # identity line, title (with correlation), theme
        geom_abline(intercept=0, slope=1, linetype="longdash", size=0.2) +
        labs(title = title, subtitle = paste("correlation:",format(r,digits=3))) + xlab(col1) + ylab(col2) +
        ## theme
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())
      
      print(p)
      
    }
  }
  return(list(r=r))
}
