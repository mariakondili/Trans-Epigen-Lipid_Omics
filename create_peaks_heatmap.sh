#!/bin/bash

#> Create heatmap from ChIPseq data of IRF5
#Help and examples : http://deeptools.readthedocs.io/en/latest/content/tools/plotHeatmap.html 

homedir="/media/kondmar/my2TB/PROJECTS/NV/IRF5_ChIPseq/"
bw_dir=$homedir"Bigwigs/"

#SAMPLES
bw_0m=$bw_dir"IRF5_WT_0m.bw"
bw_120m=$bw_dir"IRF5_WT_120m.bw"
bw_30m=$bw_dir"IRF5_WT_30m.bw"

# Tariq peaks
bed_dir=$homedir"MACS2_PeakCalling/EMTAB_2661_PE/"
Tq_120m=$bed_dir"Tariq_original_files/IRF5_WT_0h_merge.bed"
Tq_0m=$bed_dir"Tariq_original_files/IRF5_WT_2hr_merge.bed"
Tq_30m="/media/kondmar/my2TB/PROJECTS/NV/IRF5_ChIPseq/MACS2_PeakCalling/EMTAB_2033_SE/Narrow_Peaks/IRF5_30m_qval01/IRF5_30m_qval01_peaks.narrowPeak"

## with my Broad peaks of frdr=0.2: 
broad_0m=$bed_dir"BroadPeaks/cut_off_0.2/IRF5_WT_000m/Galaxy220_IRF5_WT_000m_broadPeaks.bed"
broad_30m=$homedir"MACS2_PeakCalling/EMTAB_2033_SE/Broad_cutoff_0.2/Galaxy11_IRF5_WT_030m_broadPeaks.bed"
broad_120m=$bed_dir"BroadPeaks/cut_off_0.2/IRF5_WT_120m/Galaxy221_IRF5_WT_120m_broadPeaks.bed"

# -----------------------------------------  Main command -------------------------------------# 

# Basic run No1 : computeMatrix reference-point  -S $bw_0m   -R $broad_0m   --referencePoint center -a 2000 -b 2000 -out   matrix_Peaks_0m.tab.gz

computeMatrix reference-point  -S $bw_0m   -R $broad_0m   -p max  --missingDataAsZero --referencePoint center -a 2000 -b 2000 -out   matrix_Pks_max_0m.tab.gz
computeMatrix reference-point  -S $bw_30m  -R $broad_30m  -p max  --missingDataAsZero --referencePoint center -a 2000 -b 2000 -out   matrix_Pks_max_30m.tab.gz
computeMatrix reference-point  -S $bw_120m -R $broad_120m -p max  --missingDataAsZero --referencePoint center -a 2000 -b 2000 -out   matrix_Pks_max_120m.tab.gz

deept_dir=$homedir'Deeptools_TF_heatmaps/'

plotHeatmap -m matrix_Pks_max_0m.tab.gz   -out $deept_dir'irf5_heatmap_Peaks_0m.png'   --refPointLabel  Summit --regionsLabel peaks --plotTitle 'IRF5 ChIPseq signal'
plotHeatmap -m matrix_Pks_max_30m.tab.gz  -out $deept_dir'irf5_heatmap_Peaks_30m.png'  --refPointLabel  Summit --regionsLabel peaks --plotTitle 'IRF5 ChIPseq signal'
plotHeatmap -m matrix_Pks_max_120m.tab.gz -out $deept_dir'irf5_heatmap_Peaks_120m.png' --refPointLabel  Summit --regionsLabel peaks --plotTitle 'IRF5 ChIPseq signal'

## with K-MEANS 
plotHeatmap -m matrix_Pks_max_0m.tab.gz  -out $deept_dir'kmeans_heatmaps/irf5_heatmap_Peaks_0m.png' --kmeans 5 --refPointLabel Summit  -x Distance   --plotTitle 'IRF5 ChIPseq signal'
plotHeatmap -m matrix_Pks_max_30m.tab.gz  -out $deept_dir'kmeans_heatmaps/irf5_heatmap_Peaks_30m.png' --kmeans 5 --refPointLabel Summit  -x Distance  --plotTitle 'IRF5 ChIPseq signal'
plotHeatmap -m matrix_Pks_max_120m.tab.gz -out $deept_dir'kmeans_heatmaps/irf5_heatmap_Peaks_120m.png' --kmeans 5 --refPointLabel Summit  -x Distance --plotTitle 'IRF5 ChIPseq signal'


## Run all bw.samples with Center = Ref point ,for one peak file each time : 
computeMatrix reference-point  -S $bw_0m $bw_30m  $bw_120m -R $broad_0m -p max  --missingDataAsZero --referencePoint center -a 2000 -b 2000 -out  RefPoint_3samples_pk_0m.tab.gz
plotHeatmap -m RefPoint_3samples_pk_0m.tab.gz  -out $deept_dir'htmap_Allsamples_RefPoint_0pks.png' --refPointLabel start --regionsLabel All Peaks

computeMatrix reference-point  -S $bw_0m $bw_30m  $bw_120m -R $broad_30m -p max  --missingDataAsZero --referencePoint center -a 2000 -b 2000 -out  RefPoint_3samples_pk_0m.tab.gz
plotHeatmap -m RefPoint_3samples_pk_0m.tab.gz  -out $deept_dir'htmap_Allsamples_RefPoint_30mPks.png' --refPointLabel start --regionsLabel Peaks_30m  

computeMatrix reference-point  -S $bw_0m $bw_30m  $bw_120m -R $broad_120m -p max  --missingDataAsZero --referencePoint center -a 2000 -b 2000 -out  RefPoint_3samples_pk_0m.tab.gz
plotHeatmap -m RefPoint_3samples_pk_0m.tab.gz  -out $deept_dir'htmap_Allsamples_RefPoint_120mPks.png' --refPointLabel start --regionsLabel Peaks_120m  


## Map a pool of Annot.genes (Uropa-broad) over 3 samples :

cat Transcr_coord_0m.bed Transcr_coord_30m.bed Transcr_coord_120m.bed  > pool_genes_coord_over_peaks.bed

computeMatrix reference-point  -S $bw_0m $bw_30m  $bw_120m -R pool_genes_coord_over_peaks.bed -p max  --missingDataAsZero --referencePoint center -a 2000 -b 2000 -out Matrix_3sampl_GenesCoord.tab.gz

plotHeatmap -m Matrix_3sampl_GenesCoord.tab.gz -out Heatmap_3sampl_GenesCoord.png --refPointLabel TSS --regionsLabel Annot.Genes


## Map a pool of All Broad Peaks over 3 samples
cat $broad_0m   $broad_30m $broad_120m > pool_AllPeaks.bed

computeMatrix reference-point  -S $bw_0m $bw_30m  $bw_120m -R pool_AllPeaks.bed  -p max  --missingDataAsZero --referencePoint center -a 2000 -b 2000 -out Matrix_3sampl_PoolPeaks.tab.gz
plotHeatmap -m Matrix_3sampl_PoolPeaks.tab.gz -out Heatmap_3sampl_PoolPeaks.png --refPointLabel TSS --regionsLabel Broad Peaks Pool


