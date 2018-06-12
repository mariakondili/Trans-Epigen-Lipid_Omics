#!/usr/bin/R

## Author : Maria Kondili 

## Use the Total concentration of each Lipid species (Sum over all the Carbon Lengths the lipid is measured) in order to 
##  see the overall concentration of the lipids in one treatment, comparing the two conditions WT,KO 


## Input : Average Concentration tables of each lipid calculated in script <calc_average_per_lipid.R>


#Example table of one Lipid, same for others: 
#Table : avg_cer_ccl4 
#C_len	  WT	  KO
#X16_0   3.1583 4.150
#X16_0dh 0.1083 0.1816
#X18_0   0.5033 0.7850
#X24_0   0.5200 0.9583
#X24_0dh 0.0183 0.0216
#X24_1   1.4816 2.6150
#Total   5.7900 8.7100


source("calc_average_per_lipid.R")

avg_dir <- paste0(homedir, "CCL4_AverageConcentration/Average_Tables")
avg_ccl4_tabs <- list.files(path=avg_dir,pattern="*.txt", full.names=TRUE)

lipid_names <- unlist(lapply(avg_ccl4_tabs ,function(tab_name) { 
                      avg_nm<- strsplit( basename (tab_name), ".txt")[[1]]; 
                      strsplit(avg_nm, "average_")[[1]][2]}))

avgs_ccl4 <- lapply (avg_ccl4_tabs, function(t) { 
                    avg <- read.table(t, sep="\t",as.is=TRUE, header=TRUE,row.names = 1)})

names(avgs_ccl4)<- lipids_names

totals_ccl4 <- t(sapply(lipid_names, function(L){ cbind(avgs_ccl4[[L]]["Total","WT"], avgs_ccl4[[L]]["Total","KO"])}))
colnames(totals_ccl4) <- c("WT","KO" )


write.table(totals_ccl4, file="CCL4_AverageConcentration/All_Lipids_TotalConcentr_CCl4.txt", 
            col.names = TRUE, row.names = TRUE,sep="\t", quote=FALSE)


##Do the melt in the plot.function 
#   melt_avgs_ccl4  <- melt(all_avgs_ccl4); colnames(melt_avgs_ccl4) <- c("Lipid", "Condition", "Concentration")

sem_array<-function(v) {
  sem <-sd(v)/sqrt(length(v))
  return(sem)
}


# wt and ko tables are created in script : calc_average_per_lipid.R

##> Table <wt_cer> : 
#        X16_0 X16_0dh X18_0 X24_0 X24_0dh X24_1 Total
# smpl1  2.66    0.08  0.39  0.48    0.02  1.33  4.95
# smpl2  0.39    0.00  0.02  0.11    0.02  0.10  0.64
# smpl3  4.37    0.10  0.81  0.64    0.03  2.03  7.98
# smpl4  6.59    0.21  1.27  1.19    0.03  3.62 12.91
# smpl5  0.53    0.01  0.00  0.06    0.00  0.12  0.72
# smpl6  4.41    0.25  0.53  0.64    0.01  1.69  7.54


all_wt <- list(wt_cer,wt_lpe,wt_lps,wt_ape,wt_alpe,wt_pa,wt_pe,wt_pg, wt_pi, wt_ps)
all_ko <- list(ko_cer,ko_lpe,ko_lps,ko_ape,ko_alpe,ko_pa,ko_pe,ko_pg, ko_pi, ko_ps)

sem_totals_ccl4 <- rbind ( cbind (sapply(all_wt, function(x) sem_array(x$Total),simplify=TRUE), 
                                  sapply(all_ko, function(x) sem_array(x$Total),simplify=TRUE)))
sem_totals_ccl4 <- as.data.frame(sem_totals_ccl4 ); colnames(sem_totals_ccl4 ) <- c("WT","KO"); 
rownames(sem_totals_ccl4) <- c("CER", "LPE","LPS", "aPE", "aLPE", "PA", "PE", "PG", "PI","PS")

write.table(sem_totals_ccl4, file="CCL4_AverageConcentration/SEM_AllLipids_TotalConcentr_CCl4.txt", 
            col.names = TRUE, row.names = TRUE,sep="\t", quote=FALSE)


plot.concentr.sem <- function( av_tab, sem_tab, lipid_name,treatment) { ## std_tab can be used inst.of sem_tab
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
  
  ##> For verifying if Difference of WT vs KO is significant : 
  # t <- t.test(x=aver_melt$Concentr[aver_melt$Condition =="WT"],
  #             y=aver_melt$Concentr[aver_melt$Condition =="KO"])
  # pval<- t$p.value
  
  ggplot(aver_melt, aes(x=LipidsID, y=Concentration,fill=Condition)) + 
        geom_bar(stat="identity", position="dodge") +  
        geom_errorbar(sem_lim, position="dodge")+
        labs(y="Concentration(ng/5*10^6 cells)", x="Lipids") +
        ggtitle(paste0("Sum of C-lengths Concentration of ",lipid_name, "in", treatment)) +
        scale_fill_discrete("Conditions",labels = colnames(av_tab))+
        scale_fill_brewer(palette="Paired")+
        theme(legend.position ="top") 
}


plot.concentr.sem(totals_ccl4, sem_totals_ccl4, "All Lipids","CCl4")




