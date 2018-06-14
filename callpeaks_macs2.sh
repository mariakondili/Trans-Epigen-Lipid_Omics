#!/usr/bin/bash

## Author : Maria Kondili

## Subject : Call Peaks using MACS2 in broad and Narrow mode, inputting my bam files after  trimming adapters and filtering by mapping Quality.



homedir="/home/kondmar/ChIPseq/bam/Trimmed/Filtered/"
cd $homedir

samples=(`ls *.bam`) ## ls+ path relative to working dir.
#inputsdir=(`ls bam/input/*.bam`)
ctrl=(`ls KO/*.bam`) 

outdir_narrow=$homedir'MACS2_PeakCalling/Narrow/'
outdir_broad=$homedir'MACS2_PeakCalling/Broad/'

mkdir $outdir_narrow
mkdir $outdir_broad


## PAIRED-END Reads
for i in ${!samples[@]}; do  # ${!dir[@]} : gives the indices of the elem in dir
	bam=${samples[i]}
	ko=${ctrl[i]}
	name=${bam:0:10} # capture 1st character and 10 more from bam filename

	## Narrow Peaks 
	macs2 callpeak -t $bam -c  $ko --format=BAMPE --bdg --qvalue=0.1 --mfold 10 30 --gsize 1.87e9 --name $name --outdir $outdir_narrow$name 2> $outdir_narrow$name".log"
	#OR : Broad Peaks 
	macs2 callpeak -t $bam -c  $ko --format=BAMPE --bdg --broad --broad-cutoff 0.1 --mfold 10 30 --gsize 1.87e9 --name $name --outdir $outdir_broad$name 2> $outdir_broad$name".log"

done


## Thanks to Tariq Khoyratty for suggestions on the flags to use. 