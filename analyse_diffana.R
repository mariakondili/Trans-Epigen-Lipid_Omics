#!/usr/bin/R

## Author: Maria Kondili

## Subject : Read RNA seq differential expression output and 

workdir <-"/home/kondmar/PROJECTS/FA/RNAseq/analyse_diffana/"
setwd(workdir)

tab_dir <- "/home/kondmar/PROJECTS/FA/RNAseq/DESEQ2_Output/"
genes_diffana <-read.table(paste0(tab_dir,"annotatedgenes_diffana_KO-WT_6cols.csv"), sep="\t", as.is=TRUE, header=TRUE )

## Similarly for Transcripts: 
transcr_diffana <- read.table(paste0(tab_dir,"annotatedtransc_diffana_KO-WT_6cols.csv"), sep="\t", as.is=TRUE, header=TRUE )
# !! padj column was read as "character" from csv!
transcr_diffana$padj <- as.numeric(transcr_diffana$padj)



##### V O L C A N O   P L O T ####

library(ggplot2)

cof_log10p <- 1.3
cof_log2fc <- 1.5

ggplot(transcr_diffana)+ geom_point(aes(x=log2FoldChange_KO.WT, y=-log10(pval), 
                        colour = ifelse(log2FoldChange_KO.WT > cof_log2fc & -log10(pval) > cof_log10p, "Upregulated",
                                  ifelse(log2FoldChange_KO.WT < -cof_log2fc & -log10(pval) > cof_log10p, "Downregulated", 
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


##> If gene names from ChIPseq annot are Capital Letters, 
## transform also the rna.seq gene lists to Upper Letters  >>   toupper( highly_regulated_transcr$up ) 
