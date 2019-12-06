#' @title Draw Overlap Plot
#' @description ne list overlap enrichment test + volcano plot
#' @param df something
#' @author April Kim
#' @export


### --------------------------------------------------
### gene list overlap enrichment test + volcano plot
plotOverlap <- function(df,bait,listName,listDf,drawLabel){
  
  # all proteins detected in experiment + in gene list (REMOVE BAIT)
  #preys <- unique(intersect(df$gene, listDf$gene))
  #preys <- unique(df$gene[df$gene %in% listDf$gene & df$gene != bait])
  preys <- unique(df$gene[df$gene != bait])
  
  # significant preys
  sigPreys <- unique(df$gene[df$significant & df$gene %in% preys])		
  # significant genes in gene list, limited to detected preys
  sigGenes <- listDf$gene[listDf$significant & listDf$gene %in% preys]
  
  # overlap of significant preys and significant genes
  overlap <- intersect(sigPreys,sigGenes)
  # sig preys but not sig genes
  sigPreysOnly <- setdiff(sigPreys,sigGenes)
  # sig genes but not sig preys
  sigGenesOnly <- setdiff(sigGenes,sigPreys)
  # neither sig preys nor sig genes
  neither <- setdiff(setdiff(preys,sigPreys),sigGenes)
  
  # Fisher's exact test (one-tailed)
  fisherP <- fisher.test(matrix(c(length(overlap),length(sigPreysOnly),
                                  length(sigGenesOnly),length(neither)),nrow=2),alternative="greater")$p
  
  
  # start volcano plot
  p <- ggplot(df, aes(x=logFC, y=-log10(pvalue))) +
    geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="grey") +
    xlab(bquote(log[2]*"(Fold change)")) + ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +
    
    # plot all proteins (green = significant, blue = not significant)
    geom_point(alpha=0.5, size=1.5, color=ifelse(df$significant, "springgreen3", "royalblue2")) +
    
    # label bait (red = signficant, orange = not significant)
    geom_point(subset(df, gene==bait & significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="red") + 
    geom_point(subset(df, gene==bait & !significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="orange") +
    
    # label sig genes in gene list (yellow = significant, white = not significant)
    geom_point(subset(df, gene %in% sigGenes & significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="yellow") +
    geom_point(subset(df, gene %in% sigGenes & !significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="white") +
    
    geom_point(subset(df, gene==bait | gene %in% sigGenes), mapping=aes(x=logFC, y=-log10(pvalue)),
               size=2, color="black", shape=1) +
    
    # title (with FisherP) and theme
    ggtitle(paste(listName,": ",length(sigGenes)," detected, ",length(overlap)," significant, p = ",
                  format(fisherP,digits=3),sep="")) +
    theme_bw() + theme(axis.line=element_line(color="grey"), plot.title=element_text(size=10))
  
  
  ### only draw text labels for gene list genes if drawLabel==TRUE
  if (drawLabel==TRUE) {
    p <- p + geom_text_repel(subset(df, gene==bait | gene %in% sigGenes), mapping=aes(label=gene),
                             arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.15, "lines"),
                             point.padding=unit(0.2, "lines"), color="black", size=3)
  } else {
    p <- p + geom_text_repel(subset(df, gene==bait), mapping=aes(label=gene),
                             arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.15, "lines"),
                             point.padding=unit(0.2, "lines"), color="black", size=3)
  }
  
  print(p)
}

