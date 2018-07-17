
create_volcanoPlot <- function(diffana_tab, cutoff_log10pval, cutoff_log2FC, case = 3 ) {
  
  ## Input table: A norm.counts table from diff.analysis (e.g DESeq2 : diffana) that 
  ## has been annotated with transcripts or genes and containing a column <log2FoldChange_KO.WT> and <pval>" 
  ## For having a plot with a positive y axis p-value has been used as a scale of -log10,which transforms it in positive integers.
   
  library(ggplot2)
  library(gridExtra)
  
  if (case==1) {
      # F1 : genes coloured by 1 condition --> by p-value
      condition <- with(diffana_tab, ifelse(-log10(pval) > cutoff_log10pval,
                                                       "Significant","Not Significant"))
      colour_values <- c("Significant" = "slateblue1",
                            "Not Significant" = "black")
      
  } else if (case == 2 ) {
  #F2 : genes coloured by 2 conditions --> by p-value (Signif/ Not Signif ) and the Significant ones by  log2FC
  condition <-  with(diffana_tab,  ifelse( log2FoldChange_KO.WT > cutoff_log2FC &  -log10(pval) > cutoff_log10pval, 
                                            "Signif & HighlyExpr",
                                    ifelse(log2FoldChange_KO.WT < -cutoff_log2FC  & -log10(pval) > cutoff_log10pval, 
                                            "Signif & HighlyExpr",
                                    "NotSignif & LowExpr")))
  
  colour_values <- c("Signif & HighlyExpr" = "slateblue1", 
                     "Signif & HighlyExpr" = "slateblue1",
                      "NotSignif & LowExpr" = "black"   )
  
  
  } else if (case ==3 ) {
    
  #F3 : genes coloured by 3 conditions -->  by p-value(Signif/ Not Signif)and the Significant ones : by Positive  and Negative log2FC 
  condition <- with(diffana_tab, ifelse(log2FoldChange_KO.WT > cutoff_log2FC &  
                                                         -log10(pval) > cutoff_log10pval, "Upregulated", 
                                                       ifelse(log2FoldChange_KO.WT < -cutoff_log2FC & 
                                                          -log10(pval) > cutoff_log10pval, "Downregulated", 
                                                        "Not Significant")))
  
  colour_values <- c("Upregulated"  = "darkgreen",
                      "Downregulated"   = "red",
                      "Not Significant" = "black" )
  }
  else {
    cat ( "\nERROR : Wrong case value entered ! Please indicate 1, 2 or 3.\n")
    return()
  }
  
  ### MAIN HEATMAP COMMAND ####
  ggplot(diffana_tab) + geom_point(aes(x=log2FoldChange_KO.WT, y=-log10(pval),
                                      colour =  condition)) +
    scale_color_manual(name="", values=colour_values) +
    xlim(c(-10, 10)) + ylim(c(0, 15)) + theme(legend.position = "top") + 
    xlab("log2(Fold.Change)") + ylab("-log10(p-value)")
  
}


cutoff_log10pval = 1.3
cutoff_log2FC = 1.5

## Call function 
create_volcanoPlot(transcripts_diffana, cutoff_log10pval, cutoff_log2FC,case = 3)

