#!/usr/bin/R 

## Author : Maria Kondili


filter_by_pval_log2fc <- function(diffana_tab, outdir, elem="transcript", pval_co, log2fc_co) { 
  
  ## Accept a matrix of diff.expr normalized counts with log2FC_KO.WT column calculated for each transcript/gene, by a DE-algorithm( e.g DESEQ2, edgeR)

  # Filter by pval
  pval_idx <- which(diffana_tab$pval <= pval_co)
  pval_tab <- diffana_tab[pval_idx, ]
  # Filter by log2fc
  transcr_pval_upregul  <- which(pval_tab$log2FoldChange_KO.WT >= log2fc_co)
  transcr_pval_downregul  <- which(pval_tab$log2FoldChange_KO.WT <= -log2fc_co)


  ## Keep the MGI.symbol of Genes respecting the cutoff values
  upregul_pval_names <- unique(diffana_tab$MGI.symbol[transcr_pval_upregul])
  downregul_pval_names <- unique(diffana_tab$MGI.symbol[transcr_pval_downregul])
  out_filt_names <- list('up'= upregul_pval_names, 'down'=downregul_pval_names)
  
  return(out_filt_names)
  
}


## Call the function 

transcr_names_filt <- filter_by_pval_log2fc(transcr_diffana, outdir, pval_co=0.05, log2fc_co=1.5)

cat("\nUpregulated transcripts significantly differentiated: ", length(transcr_names_filt$up))   

cat("\nDownregulated transcripts significantly differentiated: ", length(transcr_names_filt$down))


