#!/usr/bin/R

workdir <-"/media/kondmar/my2TB/PROJECTS/FA/FIBROLIV_2015/RNAseq/R_analyse_diffana/"
setwd(workdir)

# Make files executable and then call source (chmod +x )
source(paste0(workdir,"Scripts/extract_UpDownRegulElements.r"))
source(paste0(workdir,"Scripts/highly_regulated.R"))


## Call a function for genes and transcripts :
#source("R_analyse_diffana/Scripts/extract_UpDownRegulElements.r")

tab_dir <- "/home/kondmar/Documents/FA/FibroLiv2015/RNAseq/diffana/"
genes_diffana <-read.table(paste0(tab_dir,"annotatedgenes_diffana_FibroLivA2015_KO-WT_6cols.csv"), sep="\t", as.is=TRUE, header=TRUE )

## Similarly for Transcripts: 
transcr_diffana <- read.table(paste0(tab_dir,"annotatedtransc_diffana_FibroLivA2015_KO-WT_6cols.csv"), sep="\t", as.is=TRUE, header=TRUE )
# !! padj column was read as "character" from csv!
transcr_diffana$padj <- as.numeric(transcr_diffana$padj)

cutoff_FDR <- 0.1

##### V O L C A N O   P L O T ####

library(ggplot2)

ggplot(transcr_diffana)+ geom_point(aes(x=log2FoldChange_KO.WT, y=-log10(pval), 
                        colour = ifelse(log2FoldChange_KO.WT > 1.5 & -log10(pval) > 1.3, "Upregulated",
                                  ifelse(log2FoldChange_KO.WT < -1.5 & -log10(pval) > 1.3, "Downregulated", 
                                                                               "Not Significant")))) + 
                    	scale_color_manual(name="Genes selected by Cut-offs" , 
                                       	   values=c("Upregulated" = "darkgreen", 
                                                    "Downregulated" = "red", 
                                                    "Not Significant" = "black"))     +
                    	xlim(c(-10, 10)) + ylim(c(0, 10)) +
                    	xlab("log2(Fold.Change)") + ylab("-log10(p-value)") 

## Note : -log10(padj) can also be used instead of 'pval'.



## Keep the list of ALL unique Transcripts and Genes
write.table(unique(transcr_diffana$MGI.symbol), file="Input_Lists/Background_alltranscr.txt",append=FALSE,quote=FALSE,sep="\n",eol="\n", 
             col.names=FALSE,row.names = FALSE)
write.table(unique(genes_diffana$MGI.symbol), file="Input_Lists/Background_allgenes.txt",append=FALSE,quote=FALSE, sep="\n",eol="\n", 
             col.names=FALSE,row.names = FALSE)

## GO Pathways : 
#Introduce the gene lists to : http://go.princeton.edu/cgi-bin/GOTermFinder 


##> If gene names from ChIPseq annot are Capital Letters, 
## transform also the rna.seq gene lists to Upper Letters  >>   toupper( highly_regulated_transcr$up ) 

##> Find ranges of values excluding Infinity  values
# inf_val <- which(diffana_tab$log2FoldChange_KO.WT == Inf)
# range (diffana_tab$log2FoldChange_KO.WT[-inf_val])
# >> -Inf 8.032346

# negInf_val <- which(diffana_tab$log2FoldChange_KO.WT == -Inf)
# range (diffana_tab$log2FoldChange_KO.WT[-negInf_val])
# >> -6.08408      Inf


  
