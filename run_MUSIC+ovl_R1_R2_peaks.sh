#!/bin/bash

### Apply MUSIC for peak detection per replicate
### Tutorial  : https://github.com/gersteinlab/MUSIC

## Executable in :'~/MyTools/MUSIC-master/bin/' --> exported $PATH

### Run MUSIC for Enriched Regions detection on IRF5-trimmed-q10-rmDup reads : 
homedir="/media/kondmar/PROJECTS/ChIPseq/bam/Trimmed/Filtered/" 
output_dir="/media/kondmar/PROJECTS/ChIPseq/MUSIC_PeakCalling/"

WT=(ls  $homedir'WT/*.bam')
KO=(ls  $homedir'KO/*.bam')  ## Input can be used equally


## PRE-PROCESS 

for i in ${!WT[@]}; do 
	wt=${WT[i]}
	ko=${KO[i]}
	name_wt=${wt: : }
	name_ko=
	
	mkdir $output_dir$name
	samtools view $wt | MUSIC -preprocess SAM stdin $output_dir$name_wt
	samtools view $ko | MUSIC -preprocess SAM stdin $output_dir$name_ko

	mkdir $output_dir'sorted/'$name 
	mkdir $output_dir'dedup/'$name

	MUSIC -sort_reads $output_dir$name_wt  $output_dir'sorted/'$name_wt

	MUSIC -remove_duplicates  $output_dir'sorted/'$name_wt 2  $output_dir'dedup/'$name_wt
	MUSIC -remove_duplicates  $output_dir'sorted/'$name_ko 2  $output_dir'dedup/'$name_ko

	mkdir 'Peaks_chip/'$name
	MUSIC -get_TF_peaks -chip $output_dir'dedup/'$name_wt -control  $output_dir'dedup/'$name_ko -begin_l 1000 -end_l 16000 -step 1.5 -l_frag 300 -q_val 0.1 
done 

## Concatenate all chromosomes per Sample : 
mrg_peaks=$output_dir"Merged_Peaks/"
mkdir $output_dir"Merged_Peaks/"

cat 'Peaks_'$name'/*.bed' > 'All_SSERS'$name'.bed'

mrg_0_R1=$mrg_peaks'All_SSERS_WT0m_R1.bed'
mrg_0_R2=$mrg_peaks'All_SSERS_WT0m_R2.bed'
mrg_120_R1=$mrg_peaks'All_SSERS_WT120m_R1.bed'
mrg_120_R2=$mrg_peaks'All_SSERS_WT120m_R2.bed'

### INTERSECT R1+R2 peaks( unique one in one Replicate or overlaps ) 
intersectBed -wa -wb -a $mrg_0_R1 -b $mrg_0_R2 | awk 'BEGIN {OFS="\t"} {if ($2 < $11) start=$2; else start=$11} {if ($3 > $12) end=$3; else end=$12} {print  $1,start,end,($5+$14)/2}' -  > $mrg_peaks"WT_0m_R1_ovl_R2.bed"

intersectBed -wa -wb -a $mrg_120_R1 -b $mrg_120_R2 | awk 'BEGIN {OFS="\t"} {if ($2 < $11) start=$2; else start=$11} {if ($3 > $12) end=$3; else end=$12} {print  $1,start,end,($5+$14)/2}' -  > $mrg_peaks"WT_120m_R1_ovl_R2.bed"

## Find which peaks are not included in the overlaps to add them over

bedtools intersect -a  $mrg_peaks'All_SSERS_WT0m_R1.bed'   -b $mrg_peaks'WT_0m_R1_ovl_R2.bed'   -v > $mrg_peaks'WT_0m_R1_uniq.bed'
bedtools intersect -a  $mrg_peaks'All_SSERS_WT120m_R1.bed' -b $mrg_peaks'WT_120m_R1_ovl_R2.bed' -v > $mrg_peaks'WT_120m_R1_uniq.bed'

bedtools intersect -a  $mrg_peaks'All_SSERS_WT0m_R2.bed'   -b $mrg_peaks'WT_0m_R1_ovl_R2.bed'   -v > $mrg_peaks'WT_0m_R2_uniq.bed'
bedtools intersect -a  $mrg_peaks'All_SSERS_WT120m_R2.bed' -b $mrg_peaks'WT_120m_R1_ovl_R2.bed' -v > $mrg_peaks'WT_120m_R2_uniq.bed'


cat $mrg_peaks'WT_0m_R1_uniq.bed' $mrg_peaks'WT_0m_R2_uniq.bed'  $mrg_peaks"WT_0m_R1_ovl_R2.bed"   >  $mrg_peaks'WT_0m_merged+R1R2uniq.bed'

cat $mrg_peaks'WT_120m_R1_uniq.bed' $mrg_peaks'WT_0m_R2_uniq.bed'  $mrg_peaks"WT_120m_R1_ovl_R2.bed" >   $mrg_peaks'WT_120m_merged+R1R2uniq.bed'


