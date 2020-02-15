
plotVolcanoMinimal <- function(df, bait, title = '', sub1 = 'proteins detected.', sub2 ='significant.',
                               size_point = 3, size_text=3, color_alpha=0.8){
  
  
  require(ggplot2)
  require(ggrepel)
  
  nTotal <- dim(df)[1]
  nSig <- sum(df$significant==TRUE)
  
  ## inital check
  if (!all(c('logFC', 'pvalue', 'significant') %in% colnames(df))) stop('data.frame does not contain some of logFC, pvalue and signifcant.')
  
  # start volcano plot
  p <- ggplot(df, aes(x = logFC, y = -log10(pvalue))) +
    # plot all proteins (green = significant, blue = not significant)
    geom_point(alpha=0.95, size=size_point+0.3, color=ifelse(df$significant, "black", "grey"), shape=ifelse(df$significant, 1, 1)) +
    geom_point(alpha=color_alpha, size=size_point, color=ifelse(df$significant, "springgreen3", "grey")) +

    
    # label bait (red = signficant, orange = not significant)
    geom_point(subset(df, gene %in% bait & significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=size_point, color="red") + 
    geom_point(subset(df, gene %in% bait & !significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=size_point, color="orange") +
    geom_point(subset(df, gene %in% bait), mapping=aes(x=logFC, y=-log10(pvalue)), size=size_point, color="black", shape=1) +	
    geom_text_repel(subset(df, gene==bait), mapping=aes(label=gene), arrow=arrow(length=unit(0.1, 'npc')),
                    box.padding=unit(0.15, "lines"), point.padding=unit(0.2, "lines"), color="black", size=size_text) +
    
    # title (with significant count) and theme
    labs(title = title, subtitle = paste(nTotal,sub1,nSig,sub2)) + 
    geom_hline(yintercept=0, color="black") + geom_vline(xintercept=0, color="black") +
    xlab(bquote(log[2]*"[fold change]")) + ylab(bquote(-log[10]*"["*italic(.("P"))*"-value]")) +
    #scale_y_continuous(labels = function(breaks) { breaks + 1}) + 
    
    ## theme
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
  

  
  print(p)
}


