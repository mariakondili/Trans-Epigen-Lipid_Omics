#!/bin/bash

#> Create heatmap of peaks from ChIPseq data to show enriched regions around TSS/Summit 

#Help and examples : http://deeptools.readthedocs.io/en/latest/content/tools/plotHeatmap.html 

homedir="/home/kondmar/PROJECTS/ChIPseq/"
bw_dir=$homedir"Bigwigs/"

#SAMPLES
bw_0m=$bw_dir"WT_0m.bw"
bw_30m=$bw_dir"WT_30m.bw"
bw_120m=$bw_dir"WT_120m.bw"


## Broad peaks 
broad_0m=$bed_dir"BroadPeaks/WT_000m/WT_0m_broadPeaks.bed"
broad_30m=$bed_dir"BroadPeaks/WT_30m_broadPeaks.bed"
broad_120m=$bed_dir"BroadPeaks/WT_120m_broadPeaks.bed"


# -----------------------------------------  Main command -------------------------------------# 



computeMatrix reference-point  -S $bw_0m   -R $broad_0m   -p max  --missingDataAsZero --referencePoint center -a 2000 -b 2000 -out   matrix_Pks_max_0m.tab.gz
computeMatrix reference-point  -S $bw_30m  -R $broad_30m  -p max  --missingDataAsZero --referencePoint center -a 2000 -b 2000 -out   matrix_Pks_max_30m.tab.gz
computeMatrix reference-point  -S $bw_120m -R $broad_120m -p max  --missingDataAsZero --referencePoint center -a 2000 -b 2000 -out   matrix_Pks_max_120m.tab.gz

deept_dir=$homedir'Deeptools_TF_heatmaps/'

plotHeatmap -m matrix_Pks_max_0m.tab.gz   -out $deept_dir'irf5_heatmap_Peaks_0m.png'   --refPointLabel  Summit --regionsLabel peaks --plotTitle 'IRF5 ChIPseq signal'
plotHeatmap -m matrix_Pks_max_30m.tab.gz  -out $deept_dir'irf5_heatmap_Peaks_30m.png'  --refPointLabel  Summit --regionsLabel peaks --plotTitle 'IRF5 ChIPseq signal'
plotHeatmap -m matrix_Pks_max_120m.tab.gz -out $deept_dir'irf5_heatmap_Peaks_120m.png' --refPointLabel  Summit --regionsLabel peaks --plotTitle 'IRF5 ChIPseq signal'

## with K-MEANS clustering of peaks
plotHeatmap -m matrix_Pks_max_0m.tab.gz  -out $deept_dir'kmeans/irf5_heatmap_Peaks_0m.png' --kmeans 5 --refPointLabel Summit  -x Distance   --plotTitle 'IRF5 ChIPseq signal'
plotHeatmap -m matrix_Pks_max_30m.tab.gz  -out $deept_dir'kmeans/irf5_heatmap_Peaks_30m.png' --kmeans 5 --refPointLabel Summit  -x Distance  --plotTitle 'IRF5 ChIPseq signal'
plotHeatmap -m matrix_Pks_max_120m.tab.gz -out $deept_dir'kmeans/irf5_heatmap_Peaks_120m.png' --kmeans 5 --refPointLabel Summit  -x Distance --plotTitle 'IRF5 ChIPseq signal'


##>  Results show great similarity with Summit or K-means method


