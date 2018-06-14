#!/usr/bin/R

## Author : Maria Kondili

## Subject : Remove Predicted genes(Gm*, *Rik) from List of gene names given

remove_pred_genes <- function(GeneNames_List, outdir, outfile_name){

	## Input: GeneNamesList as an 1D R-array, that has been read from the table of 
	GeneNames_List <- unique(GeneNames_List)

	GeneList_noRik <- chip_rna_ovl[which(!grepl("*Rik", GeneNames_List) )]
	##Verify :  grep("*Rik",ChIP_RNA_ovl_30m_noRik) >> integer (0)

	GeneList_noRik_noGm <- GeneList_noRik[which(!grepl("Gm*", GeneList_noRik))]

	write.table(GeneList_noRik_noGm, paste0(outdir,outfile_name,".txt"), 
            	    sep="\t", col.names = FALSE, row.names = FALSE, quote=FALSE)

	return(GeneList_noRik_noGm)
}


#### Do the filtering for a list of DE GeneNames coming from RNA-seq
rnaseqdir <- "/home/kondmar/RNAseq/diffana/Genes_log2FC_1.5/"

genes_outdir <- paste0(rnaseqdir, 'Filtered_Gm+Rik_genes/')

#Create the R-array object of the one column-GenesList
Upregul_GenesList <- read.table(paste0(rnaseqdir,"Upregulated_Genes.txt"), sep="\t", header=FALSE, as.is=TRUE)[,1]

upregGen_noGmRik <-remove_pred_genes(Upregul_GenesList, genes_outdir, "UpregGenes_NoGmRik")

##> Output: New GenesList will be written in genes_outdir.
