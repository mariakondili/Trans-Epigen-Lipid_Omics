#!/bin/bash

### Apply MUSIC for peak detection per replicate
### Tutorial  : https://github.com/gersteinlab/MUSIC

## Executable in :'~/MyTools/MUSIC-master/bin/' --> exported $PATH

### Run MUSIC for Enriched Regions detection on IRF5-trimmed-q10-rmDup reads : 
homedir="/media/kondmar/my2TB/PROJECTS/NV/IRF5_ChIPseq/bam/BAM_Trim_q10_RmDup/"
output_dir="/media/kondmar/my2TB/PROJECTS/NV/IRF5_ChIPseq/MUSIC_PeakCalling/"

cd $homedir
WT=(ls WT/*.bam)
WT30='WT/SE/WT_030m_R1_sorted_rmDup.bam'
KO=(ls KO/*.bam)
input_30='Input_SE/IRF5_Input_SE_sorted_q10_rmDup.bam'

## PRE-PROCESS 
mkdir $output_dir'chip_0m_R1'
samtools view ${WT[1]} | MUSIC -preprocess SAM stdin $output_dir'chip_0m_R1'
mkdir $output_dir'chip_0m_R2'
samtools view ${WT[2]} | MUSIC -preprocess SAM stdin $output_dir'chip_0m_R2'

samtools view ${WT[3]} | MUSIC -preprocess SAM stdin $output_dir'chip_120m_R1'
samtools view ${WT[4]} | MUSIC -preprocess SAM stdin $output_dir'chip_120m_R2'

mkdir $output_dir'chip_30m'
samtools view $WT30  | MUSIC -preprocess SAM stdin $output_dir'chip_30m'
samtools view $input_30  | MUSIC -preprocess SAM stdin $output_dir'input_30m'

samtools view ${KO[1]} | MUSIC -preprocess SAM stdin $output_dir'KO_0m_R1'
samtools view ${KO[2]} | MUSIC -preprocess SAM stdin $output_dir'KO_0m_R2'
samtools view ${KO[3]} | MUSIC -preprocess SAM stdin $output_dir'KO_120m_R1'
samtools view ${KO[4]} | MUSIC -preprocess SAM stdin $output_dir'KO_120m_R2'

mkdir sorted/chip_0m;mkdir sorted/chip_120m; mkdir sorted/chip_30m
mkdir sorted/KO_0m; mkdir sorted/KO_120m;mkdir sorted/input_30m

MUSIC -sort_reads $output_dir'chip_0m_R1'   $output_dir'sorted/chip_0m_R1'
MUSIC -sort_reads $output_dir'chip_0m_R2'   $output_dir'sorted/chip_0m_R2'
#120m
MUSIC -sort_reads $output_dir'chip_120m_R1' $output_dir'sorted/chip_120m_R1'
MUSIC -sort_reads $output_dir'chip_120m_R2' $output_dir'sorted/chip_120m_R2'

MUSIC -sort_reads $output_dir'chip_30m' $output_dir'sorted/chip_30m'
MUSIC -sort_reads $output_dir'input_30m' $output_dir'sorted/Input_30m'

MUSIC -sort_reads $output_dir'KO_0m_R1' $output_dir'sorted/KO_0m_R1'
MUSIC -sort_reads $output_dir'KO_0m_R2' $output_dir'sorted/KO_0m_R2'

MUSIC -sort_reads $output_dir'KO_120m_R1' $output_dir'sorted/KO_120m_R1'
MUSIC -sort_reads $output_dir'KO_120m_R2' $output_dir'sorted/KO_120m_R2'


MUSIC -remove_duplicates  $output_dir'sorted/chip_0m_R1' 2  $output_dir'dedup/chip_0m_R1'
MUSIC -remove_duplicates  $output_dir'sorted/chip_0m_R2' 2  $output_dir'dedup/chip_0m_R2'

MUSIC -remove_duplicates  $output_dir'sorted/chip_120m_R1' 2  $output_dir'dedup/chip_120m_R1'
MUSIC -remove_duplicates  $output_dir'sorted/chip_120m_R2' 2  $output_dir'dedup/chip_120m_R2'

MUSIC -remove_duplicates $output_dir'sorted/chip_30m' 2  $output_dir'dedup/chip_30m'
MUSIC -remove_duplicates $output_dir'sorted/Input_30m' 2 $output_dir'dedup/Input_30m'

MUSIC -remove_duplicates $output_dir'sorted/KO_0m_R1' 2 $output_dir'dedup/KO_0m_R1'
MUSIC -remove_duplicates $output_dir'sorted/KO_0m_R2' 2 $output_dir'dedup/KO_0m_R2'

MUSIC -remove_duplicates $output_dir'sorted/KO_120m_R1' 2  $output_dir'dedup/KO_120m_R1'
MUSIC -remove_duplicates $output_dir'sorted/KO_120m_R2' 2  $output_dir'dedup/KO_120m_R2'


cd /media/kondmar/my2TB/PROJECTS/NV/IRF5_ChIPseq/MUSIC_PeakCalling/dedup/

#MUSIC -get_TF_peaks -chip chip -control input -begin_l 1000 -end_l 16000 -step 1.5 -l_frag 75 	-q_val 0.1

mkdir Peaks_chip_0m_R1/
cd Peaks_chip_0m_R1/

MUSIC -get_TF_peaks -chip $output_dir'dedup/chip_0m_R1' -control  $output_dir'dedup/KO_0m_R1' -begin_l 1000 -end_l 16000 -step 1.5 -l_frag 300 -q_val 0.1 
MUSIC -get_TF_peaks -chip $output_dir'dedup/chip_0m_R2' -control  $output_dir'dedup/KO_0m_R2' -begin_l 1000 -end_l 16000 -step 1.5 -l_frag 300 -q_val 0.1

mkdir Peaks_chip_120m_R1, _R2

MUSIC -get_TF_peaks -chip $output_dir'dedup/chip_120m_R1' -control  $output_dir'dedup/KO_120m_R1' -begin_l 1000 -end_l 16000 -step 1.5 -l_frag 300 -q_val 0.1
MUSIC -get_TF_peaks -chip $output_dir'dedup/chip_120m_R2' -control  $output_dir'dedup/KO_120m_R2' -begin_l 1000 -end_l 16000 -step 1.5 -l_frag 300 -q_val 0.1

mkdir Peaks_chip_30m
cd Peaks_chip_30m
MUSIC -get_TF_peaks -chip $output_dir'dedup/chip_30m' -control  $output_dir'dedup/Input_30m' -begin_l 1000 -end_l 16000 -step 1.5 -l_frag 300 -q_val 0.1


## Concatenate all chromosomes per Sample : 
mrg_peaks=$output_dir"Merged_Peaks/"
mkdir $output_dir"Merged_Peaks/"
cat Peaks_chip_0m_R1/*.bed >  All_SSERS_WT0m_R1.bed
cat Peaks_chip_120m_R1/*.bed > All_SSERS_WT30m.bed
cat Peaks_chip_30m/*.bed >  All_SSERS_WT120m_R1.bed

mrg_0_R1=$mrg_peaks'All_SSERS_WT0m_R1.bed'
mrg_0_R2=$mrg_peaks'All_SSERS_WT0m_R2.bed'
mrg_120_R1=$mrg_peaks'All_SSERS_WT120m_R1.bed'
mrg_120_R2=$mrg_peaks'All_SSERS_WT120m_R2.bed'

### INTERSECT R1+R2 peaks( unique one in one Replicate or overlaps ) 
intersectBed -wa -wb -a $mrg_0_R1 -b $mrg_0_R2 | awk 'BEGIN {OFS="\t"} {if ($2 < $11) start=$2; else start=$11} {if ($3 > $12) end=$3; else end=$12} {print  $1,start,end,($5+$14)/2, $6,$7, $8, $9}' -  > $mrg_peaks"WT_0m_R1_ovl_R2.bed"

intersectBed -wa -wb -a $mrg_120_R1 -b $mrg_120_R2 | awk 'BEGIN {OFS="\t"} {if ($2 < $11) start=$2; else start=$11} {if ($3 > $12) end=$3; else end=$12} {print  $1,start,end,($5+$14)/2, $6,$7, $8, $9}' -  > $mrg_peaks"WT_120m_R1_ovl_R2.bed"

## Find which peaks are not included in the overlaps to add them over

bedtools intersect -a  $mrg_peaks'All_SSERS_WT0m_R1.bed'   -b $mrg_peaks'WT_0m_R1_ovl_R2.bed'   -v > $mrg_peaks'WT_0m_R1_uniq.bed'
bedtools intersect -a  $mrg_peaks'All_SSERS_WT120m_R1.bed' -b $mrg_peaks'WT_120m_R1_ovl_R2.bed' -v > $mrg_peaks'WT_120m_R1_uniq.bed'

bedtools intersect -a  $mrg_peaks'All_SSERS_WT0m_R2.bed'   -b $mrg_peaks'WT_0m_R1_ovl_R2.bed'   -v > $mrg_peaks'WT_0m_R2_uniq.bed'
bedtools intersect -a  $mrg_peaks'All_SSERS_WT120m_R2.bed' -b $mrg_peaks'WT_120m_R1_ovl_R2.bed' -v > $mrg_peaks'WT_120m_R2_uniq.bed'

cat $mrg_peaks'WT_0m_R1_uniq.bed' $mrg_peaks'WT_0m_R2_uniq.bed'  $mrg_peaks"WT_0m_R1_ovl_R2.bed"   >  $mrg_peaks'WT_0m_merged+R1R2uniq.bed'
cat $mrg_peaks'WT_120m_R1_uniq.bed' $mrg_peaks'WT_0m_R2_uniq.bed'  $mrg_peaks"WT_120m_R1_ovl_R2.bed" >   $mrg_peaks'WT_120m_merged+R1R2uniq.bed'


