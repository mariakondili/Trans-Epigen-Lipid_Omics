#!/usr/bin/R
# Author : Maria Kondili
homedir <- "home/PROJECTS/Lipidomics/Lipids_AverageConcentration/"
setwd(homedir)
library(ggplot2)
library(reshape2)


#######  Calculate Average for each lipid  species #######

# Lipids concentration must be saved in smaller tables per condition: WT,KO (columns) for all samples (lines)
#Cer	WT	KO	
#S1	.	.
#S2	.	.
#S3	.	.
#S4	.	.

##> Read tables 

ccl4_files <- paste0(homedir, "CCL4_input_tables/")
wt_cer <- read.table(paste0(ccl4_files,"WT_Cer.txt"), sep="\t", as.is=TRUE, header=TRUE,row.names = 1)
ko_cer <- read.table(paste0(ccl4_files,"KO_Cer.txt"), sep="\t", as.is=TRUE, header=TRUE,row.names = 1)

## sum up the total of each sample on Carbon-length concentr. in a global total 
total_cer <- mean(wt_cer$Total); 
sem_cer_wt <- sem_array(wt_cer$Total)


##> Calculate AVERAGE PER COLUMN AND Save in new table WT + KO 
avg_cer <- as.data.frame(cbind(colMeans(wt_cer), colMeans(ko_cer))); 
colnames(avg_cer)<- c("WT","KO")

#rownames(avg_cer)[2] <-"X16_0dh"  
#rownames(avg_cer)[5] <-"X24_0dh"


## Write down in a file : 
write.table(avg_cer, paste0(homedir,"Average_Tables/average_Cer.txt"), 
			    sep="\t", quote=FALSE, col.names = T, row.names = T)


##> Calculate STD per Column(samples)
stdev_cer <- cbind( sapply(wt_cer, sd, 2), sapply(ko_cer, sd, 2)); colnames(stdev_cer) <- c("WT","KO")


####### Calculate STD Error of the Mean  = S E M #####

calc_sem <- function(wt_tab,ko_tab) {
  ## input :2-D TABLES from lipids with Rows=samples, Columns =C-lengths
  ## Calculates Std Error of Means for a Lipid on each C-length and for each condition:WT,KO
  sem_tab <- cbind( sapply(wt_tab, sd, 2)/sqrt(nrow(wt_tab)),
                    sapply(ko_tab, sd, 2)/sqrt(nrow(ko_tab)) )
  colnames(sem_tab) <- c("WT","KO")
  ## Returns SEM for each C-length (Lines) on WT,KO =columns
  return(sem_tab)
}

# Call function
sem_cer <- calc_sem(wt_cer, ko_cer)

source("plot.concentr.sem.R")
outdir <- paste0(homedir, "AverConcentr_Plot/") 

if  (!dir.exists(outdir) ) { 
	dir.create(outdir )  
}


## Function 
plot.concentr.sem <- function( av_tab, sem_tab, lipid_name,treatment) { 
				## std_tab can be used inst.of sem_tab
  require(reshape2)            
  library(ggplot2)
  aver_melt <- melt(av_tab)
  sem_melt <- melt(sem_tab)
  
  ## columns created :
  # lipid_id -> rowname 
  # variable ->  stacked : average-wt, average-ko
  # value -> concentr.average

  colnames(aver_melt)<- c("LipidsID","Condition","Concentration" ) 

  sem_lim <- aes(ymin=as.numeric(aver_melt$Concentration) - as.numeric(sem_melt$value) ,
                 ymax=as.numeric(aver_melt$Concentration) + as.numeric(sem_melt$value))
  

  ggplot(aver_melt, aes(x=LipidsID, y=Concentration,fill=Condition)) + 
        geom_bar(stat="identity", position="dodge") +  
        geom_errorbar(sem_lim, position="dodge")+
        #theme(axis.text.x = element_text(angle = 40, hjust = 1))+
        labs(y="Concentration(ng/5*10^6 cells)", x="Lipids") +
        ggtitle(paste0("Average Concentration of ",lipid_name)) +
        scale_fill_discrete("Conditions",labels = colnames(av_tab))+
        scale_fill_brewer(palette="Paired")+
        theme(legend.position ="top")
      
}

## Call function for plotting 
plot.concentr.sem(avg_cer, sem_cer,"Cer")
  




