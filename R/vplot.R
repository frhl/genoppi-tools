
#  a function that is goin to replace plotVolcano and plotOverlap


vplot <- function(data, bait, overlap = NULL, title = '', subtitle = ''){
  
  require(ggplot2)
  require(ggrepel)
  
  nTotal <- dim(df)[1]
  nSig <- sum(df$significant==TRUE)
  
  # generate statistics for enrichement
  statistics <- enrichment(df, bait, reference)
  
  # minimal volcano plot
  plt <- ggplot(df, aes(x = logFC, y = -log10(pvalue))) +
    geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="grey") +
    xlab(bquote(log[2]*"(Fold change)")) + ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +
    # plot all proteins (green = significant, blue = not significant)
    # label bait (red = signficant, orange = not significant)
    geom_point(alpha=0.5, size=1.5, color=ifelse(df$significant, "springgreen3", "royalblue2")) +
    geom_point(subset(df, gene==bait & significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="red") + 
    geom_point(subset(df, gene==bait & !significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="orange")
  
  # draw overlapping points overlapping with genelist
  plt <- plt + geom_point(subset(df, gene %in% statistics$sigGenes & significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color='yellow') +
    geom_point(subset(df, gene %in% statistics$sigGenes & !significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="white") +
    
    geom_point(subset(df, gene==bait | gene %in% statistics$sigGenes), mapping=aes(x=logFC, y=-log10(pvalue)),
               size=2, color="black", shape=1) +
    
    #geom_point(subset(df, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="black", shape=1) +	
    #geom_text_repel(subset(df, gene==bait), mapping=aes(label=gene), arrow=arrow(length=unit(0.015, 'npc')),
    #                box.padding=unit(0.15, "lines"), point.padding=unit(0.2, "lines"), color="black", size=3) +
    
    # title (with significant count) and theme
    labs(title = title, subtitle = paste(nTotal,sub1,nSig,sub2)) +
    theme_bw() + theme(axis.line=element_line(color="grey")) + ggstamp()
  
  
  
}