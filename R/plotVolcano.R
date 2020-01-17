#' @title plot volcano
#' @description Volcano plot
#' @param df a data.frame containing columns logFC, significant and p-value.
#' @param bait what bait was used?  a string.
#' @param title should the ggplot ahve a title?
#' @param sub1 The first part of a subtitle that is to be generated.
#' @param sub2 The second part of the subtitle.
#' @author April Kim
#' @family genoppi
#' @export



plotVolcano <- function(df, bait, title = '', sub1 = 'proteins detected.', sub2 ='significant.'){
  
  require(ggplot2)
  require(ggrepel)
  
  nTotal <- dim(df)[1]
  nSig <- sum(df$significant==TRUE)
  
  ## inital check
  if (!all(c('logFC', 'pvalue', 'significant') %in% colnames(df))) stop('data.frame does not contain some of logFC, pvalue and signifcant.')
  
  # start volcano plot
  p <- ggplot(df, aes(x = logFC, y = -log10(pvalue))) +
    geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="grey") +
    xlab(bquote(log[2]*"(Fold change)")) + ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +
    
    # plot all proteins (green = significant, blue = not significant)
    geom_point(alpha=0.5, size=1.5, color=ifelse(df$significant, "springgreen3", "royalblue2")) +
    
    # label bait (red = signficant, orange = not significant)
    geom_point(subset(df, gene %in% bait & significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="red") + 
    geom_point(subset(df, gene %in% bait & !significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="orange") +
    geom_point(subset(df, gene %in% bait), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="black", shape=1) +	
    geom_text_repel(subset(df, gene==bait), mapping=aes(label=gene), arrow=arrow(length=unit(0.015, 'npc')),
                    box.padding=unit(0.15, "lines"), point.padding=unit(0.2, "lines"), color="black", size=3) +
    
    # title (with significant count) and theme
    labs(title = title, subtitle = paste(nTotal,sub1,nSig,sub2)) +
    theme_bw() + theme(axis.line=element_line(color="grey")) + ggstamp()
  
  print(p)
}