#!/usr/bin/R

## Author : Maria Kondili 

## Use the Total concentration of each Lipid species (Sum over all the Carbon Lengths the lipid is measured) in order to 
##  see the overall concentration of the lipids in one treatment, comparing the two conditions WT,KO 


## Input : From the Concentration tables of each condition,
## calculate the Sum of all carbon lengths and save it in the table as a new column.

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

totals_ccl4 <- rbind (cbind (avg_cer_ccl4["Total", "WT"], avg_cer_ccl4["Total","KO"]),
                        cbind (avg_lpe_ccl4["Total","WT"],avg_lpe_ccl4["Total","KO"]),
                        cbind (avg_lps_ccl4["Total","WT"],avg_lps_ccl4["Total", "KO"]),
                        cbind (avg_ape_ccl4["Total","WT"],avg_ape_ccl4["Total","KO"]),
                        cbind (avg_alpe_ccl4["Total","WT"],avg_alpe_ccl4["Total","KO"]),
                        cbind (avg_pa_ccl4["Total","WT"],avg_pa_ccl4["Total", "KO"]),
                        cbind (avg_pe_ccl4["Total","WT"],avg_pe_ccl4["Total","KO"]),
                        cbind (avg_pg_ccl4["Total","WT"],avg_pg_ccl4["Total","KO"]),
                        cbind (avg_pi_ccl4["Total","WT"],avg_pi_ccl4["Total", "KO"]),
                        cbind (avg_ps_ccl4["Total","WT"],avg_ps_ccl4["Total","KO"]))

rownames(totals_ccl4) <- c("CER", "LPE","LPS", "aPE", "aLPE", "PA", "PE", "PG", "PI", "PS")
colnames(totals_ccl4)<- c("WT","KO")

write.table(totals_ccl4, file="CCL4_AverageConcentration/All_Lipids_TotalConcentr_CCl4.txt", 
            col.names = TRUE, row.names = TRUE,sep="\t", quote=FALSE)


##Do the melt in the plot.function 
#   melt_avgs_ccl4  <- melt(all_avgs_ccl4); colnames(melt_avgs_ccl4) <- c("Lipid", "Condition", "Concentration")

sem_array<-function(v) {
  sem <-sd(v)/sqrt(length(v))
  return(sem)
}

sem_totals_ccl4 <- rbind(cbind(sem_array(wt_cer$Total),sem_array(ko_cer$Total)), 
                       cbind(sem_array(wt_lpe$Total), sem_array(ko_lpe$Total)),
                       cbind(sem_array(wt_lps$Total), sem_array(ko_lps$Total)),
                       cbind(sem_array(wt_ape$Total), sem_array(ko_ape$Total)),
                       cbind(sem_array(wt_alpe$Total),sem_array(ko_alpe$Total)),
                       cbind(sem_array(wt_pa$Total), sem_array(ko_pa$Total)),
                       cbind(sem_array(wt_pe$Total), sem_array(ko_pe$Total)),
                       cbind(sem_array(wt_pg$Total), sem_array(ko_pg$Total)),
                       cbind(sem_array(wt_pi$Total), sem_array(ko_pi$Total)),
                       cbind(sem_array(wt_ps$Total), sem_array(ko_ps$Total)) )



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
  
  # t <- t.test(x=aver_melt$Concentr[aver_melt$Condition =="WT"],
  #             y=aver_melt$Concentr[aver_melt$Condition =="KO"])
  # pval<- t$p.value
  
  ggplot(aver_melt, aes(x=LipidsID, y=Concentration,fill=Condition)) + 
        geom_bar(stat="identity", position="dodge") +  
        geom_errorbar(sem_lim, position="dodge")+
        labs(y="Concentration(ng/5*10^6 cells)", x="Lipids") +
        #ggtitle(paste0("Sum of C-lengths Concentration of ",lipid_name, "in", treatment)) +
        scale_fill_discrete("Conditions",labels = colnames(av_tab))+
        scale_fill_brewer(palette="Paired")+
        theme(legend.position ="top") 
}


plot.concentr.sem(totals_ccl4, sem_totals_ccl4, "All Lipids","CCl4")




