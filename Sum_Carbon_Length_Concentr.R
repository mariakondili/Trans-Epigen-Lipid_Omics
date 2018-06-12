#!/usr/bin/R
## Author : Maria Kondili

##>Subject : Design in lines the progression of Carbon-lengths for lipids with 32-42 C in order to show Elongation 
## elongation_plot.R
workdir<- "/media/kondmar/PROJECTS/Lipidomics/SUM_per_CarbonLength/"

ccl4_dir <- paste0(workdir,"CCl4/input_tabs/")

# Tables are in form : Samples(rows) X carbon_length(column)
### Example Table : 
#        C_32       C_34       C_36      C_38       C_40       C_42
#S1 0.68849496  5.9845025  29.080874 15.232200 16.4702555 0.25709567
#S2 0.00377096  0.3448007   3.059212  1.007609  0.8246590 0.05556182
#S3 1.48390446 10.4102577  45.591192 26.246138 28.4443701 0.50009281
#S4 2.06411963 19.7344609 104.194986 51.910519 51.1887148 1.72562701
#S5 0.01791238  0.4351468   3.674048  1.060974  0.9326841 0.05418961
#S6 1.32143187 12.4239220  61.397574 31.228158 30.6786562 0.68466168


long_lipids<-c("aPE", "PA","PE","PG", "PI","PS")

WT_tabs <- lapply (long_lipids, function(lip_nm) { read.table(paste0(ccl4_dir, "WT_", lip_nm, ".txt"), 
                                                              sep="\t",as.is=TRUE, na.strings = "NA",header=TRUE,row.names=1)
});names(WT_tabs)<- long_lipids


KO_tabs <- lapply (long_lipids, function(lip_nm) { read.table(paste0(ccl4_dir, "KO_", lip_nm, ".txt"), 
                                                              sep="\t",as.is=TRUE, na.strings = "NA",header=TRUE, row.names=1)
});names(KO_tabs) <- long_lipids


rbind_tab <- function(fun2bind=colMeans, tab, lipids_list) { 
  
    fun_per_lipid <- lapply(lipids_list,function(L) {fun2bind(tab[[L]])} )
    names(fun_per_lipid) <- lipids_list
    fun_per_lipid  <- as.data.frame(fun_per_lipid )
    return(fun_per_lipid)
}


###### Calculate Mean +Sem for each Lipid table per Carbon Length  ######
source(paste0(homedir,"R_functions/calc_sem.R"))

mean_perClen_WT <- rbind_tab (colMeans, WT_tabs,long_lipids )
sem_perClen_WT <- rbind_tab(fun2bind = calc_sem_1tab, WT_tabs, long_lipids  )

# KO 
mean_perClen_KO <- rbind_tab (colMeans, KO_tabs,long_lipids )
sem_perClen_KO <- rbind_tab(fun2bind = calc_sem_1tab, KO_tabs, long_lipids  )


library(reshape2)
## Create Reformatted Matrices for 
mt_Mean_Clen_WT <- melt(as.matrix(mean_perClen_WT))
mt_sem_Clen_WT <- melt(as.matrix(sem_perClen_WT))

mt_Mean_Clen_KO <- melt(as.matrix(mean_perClen_KO))
mt_sem_Clen_KO <- melt(as.matrix(sem_perClen_KO))

colnames(mt_Mean_Clen_WT ) <- colnames(mt_Mean_Clen_KO) <- c("Lipids" , "CarbonLength", "Concentration")
colnames(mt_sem_Clen_WT )  <- colnames(mt_sem_Clen_KO )  <- c("Lipids" , "CarbonLength", "sem")
 

library(ggplot2) ; library(scales) 

sem_Clen_WT_limits <- aes(ymin=as.numeric(mt_Mean_Clen_WT$Concentration) - as.numeric(mt_sem_Clen_WT$sem) ,
                          ymax=as.numeric(mt_Mean_Clen_WT$Concentration) + as.numeric(mt_sem_Clen_WT$sem) )

sem_Clen_KO_limits <- aes(ymin=as.numeric(mt_Mean_Clen_KO$Concentration) - as.numeric(mt_sem_Clen_KO$sem) ,
                          ymax=as.numeric(mt_Mean_Clen_KO$Concentration) + as.numeric(mt_sem_Clen_KO$sem) )


###  Sum WT & KO Lipids families Concentration TO show One Line as a global image.

sum_Lipids_WT <- melt(colSums(mean_perClen_WT[,-6])) #C42 will be removed bcs of NAs
sum_Lipids_KO <- melt(colSums(mean_perClen_KO[,-6]))
ColSum_perClen <- cbind(sum_Lipids_WT, sum_Lipids_KO )

colnames(ColSum_perClen) <- c("WT","KO")
ColSum_perClen[,"CarbonLength"]<- rownames(mt_Mean_Clen_WT)

mt_sumLipids <- melt(ColSum_perClen, 
			id.vars = "CarbonLength", 
			variable.name = "Conditions",
			value.name = "Concentration")

### Ready for plotting : 

library(scales)

ggplot(data=mt_sumLipids, aes(x=CarbonLength, y=Concentration, group=Conditions,color=Conditions)) + 
            geom_point(size=2) + geom_line(linetype=1,size=1.4) +  
            labs(title="Sum of Lipids Concentration per Carbon Length",y="Sum.Concentration") + 
            scale_y_continuous(breaks=pretty_breaks(n=8))+
            theme(legend.position="top")

