#' @title Draw Overlap Plot
#' @description ne list overlap enrichment test + volcano plot
#' @param df something
#' @author April Kim
#' @family genoppi
#' @export


### --------------------------------------------------
### gene list overlap enrichment test + volcano plot
plotOverlapMinimal <- function(df, bait, reference, title = '', subtitle = NULL, drawLabel = T,
                        col.genelist.sig = 'yellow', inweb_enrichment_calculate = F,
                        size_point = 3, size_text=3, color_alpha=0.8, 
                        plot.legend = T, drawLabelOnlySignificant = T,
                        mod.arrowhead = 0.5, mod.arrowlength = 3,
                        bait.cols = c('red', 'orange')){
  
  require(ggplot2)
  require(ggrepel)
  
  # If a df with only genes are inputted assume
  # that the user would like to just overlay genes
  #if (ncol(reference) == 1){ 
  #  reference <- data.frame(gene=reference$gene, significant=TRUE)
  #  warn('[plotOverlap] no "significant"-column, assuming al#l genes significant.')
  #}
  
  # generate statistics for enrichement
  if (inweb_enrichment_calculate){
    statistics <- enrichment_inweb(df, bait, reference)
  } else{
    statistics <- enrichment(df, bait, reference)
  }

  
  
  if ('draw' %nin% colnames(reference)) reference$draw <- TRUE
  
  subset1 <- merge(subset(df, gene %in% statistics$sigGenes & significant), reference)
  subset2 <- merge(subset(df, gene %in% statistics$sigGenes & !significant), reference)
  
  # start volcano plot
  p <- ggplot(df, aes(x=logFC, y=-log10(pvalue))) +
    geom_hline(yintercept=0, color="black") + geom_vline(xintercept=0, color="black") +
    xlab(bquote(log[2]*"[fold change]")) + ylab(bquote(-log[10]*"["*italic(.("P"))*"-value]")) +
    
    # plot all proteins (green = significant, blue = not significant)
    geom_point(alpha=0.95, size=size_point+0.3, color=ifelse(df$significant, "black", "grey"), shape=ifelse(df$significant, 1, 1)) +
    geom_point(alpha=color_alpha, size=size_point, color=ifelse(df$significant, "springgreen3", "grey")) +
    #geom_point(alpha=color_alpha, size=size_point, color=ifelse(df$significant, "springgreen3", "grey")) +
    
    #geom_point(subset(df, gene %in% bait), mapping=aes(x=logFC, y=-log10(pvalue)), size=size_point, color="black", shape=1) +
    
    # label sig genes in gene list (yellow = significant, white = not significant)
    geom_point(subset1, mapping=aes(x=logFC, y=-log10(pvalue), colour = dataset), size=size_point) +
    #geom_point(subset1, mapping=aes(x=logFC, y=-log10(pvalue)), size=size_point, colour = col.genelist.sig) +
    #scale_colour_manual(values=setNames(c('orange','blue'), author)) +
    #rscale_colour_identity("author", breaks=c('orange','blue'), guide="legend") +
    geom_point(subset2, mapping=aes(x=logFC, y=-log10(pvalue)), size=size_point, color="white") +
    geom_point(subset(df, gene %in% bait | gene %in% statistics$sigGenes), mapping=aes(x=logFC, y=-log10(pvalue)),
               size=size_point, color="black", shape=1) +
    
    # label bait (red = signficant, orange = not significant)
    geom_point(subset(df, gene %in% bait & significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=size_point, color=bait.cols[1]) + 
    geom_point(subset(df, gene %in% bait & !significant), mapping=aes(x=logFC, y=-log10(pvalue)), size=size_point, color=bait.cols[2]) +
    
    # title (with statistics$fisherP) and theme
    labs(title = title,
         subtitle =  ifelse(is.null(subtitle), 
                            paste(length(statistics$sigGenes)," detected. ",length(statistics$overlap),
                                  " significant. p-value = ", format(statistics$fisherP,digits=3),sep=""),
                            subtitle))+ 
    
    #ggtitle( +
    #theme_minimal()
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), legend.position=c(.15,.15),
          legend.key = element_rect(colour = "transparent", fill = "white"))
  
  ### only draw text labels for gene list genes if drawLabel==TRUE
  if (drawLabel==TRUE) {
    p <- p + geom_point(subset(df, (gene %in% bait | gene %in% statistics$sigGenes) & (gene %in% reference[reference$draw, ]$gene) | gene %in% bait), mapping=aes(x=logFC, y=-log10(pvalue)),
                            size=size_point, color="black", shape=1)
    p <- p + geom_text_repel(subset(df, (gene==bait | gene %in% statistics$sigGenes) & (gene %in% reference[reference$draw, ]$gene) | gene %in% bait), mapping=aes(label=gene),
                            arrow=arrow(length=unit(0.015*mod.arrowhead, 'npc')), box.padding=unit(0.15*mod.arrowlength, "lines"),
                            point.padding=unit(0.2, "lines"), color="black", size=size_text)
  } else {
    p <- p + geom_text_repel(subset(df, gene==bait), mapping=aes(label=gene),
                             arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.15, "lines"),
                             point.padding=unit(0.2, "lines"), color="black", size=size_text)
  }
  
  if (TRUE){
    f <- data.frame(table(reference$dataset, reference$color))
    dataset_authors <- as.vector(f$Var1[as.logical(f$Freq)])
    dataset_colors <- as.vector(f$Var2[as.logical(f$Freq)])
    p <- p + scale_color_manual(values = dataset_colors)
  }
  
  
  print(p)
}

