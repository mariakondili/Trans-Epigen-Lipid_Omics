#!/bin/bash

## Author : Maria Kondili 

## Subject :  First basic steps to do on bam files of ChIPseq experiment ,to facilitate Downstream Analysis
##            Sort and index a bam file is necessary for all NGS data after alignment. 
##            bamCoverage is done for obtaining visualisable files from which we can verify the aligment signal and observe the peaks/enriched regions of chromatin state.


inputdir='/home/kondmar/PROJECTS/ChIPseq/bam/'
sorted_dir=$inputdir'sorted/'
mkdir $sorted_dir
outdir='/home/kondmar/PROJECTS/ChIPseq/Bigwigs/'



for f in $(ls $inputdir/*.bam); do
	
	basenm=`echo $f | rev | cut -d"/" -f 1 |rev | sed 's/.bam//'`  ## rev = reads backwards, from end to beginning, so 1st field is the file.bam
	
	# Sort BAM file
	echo $basenm
	samtools sort -O BAM $f -o $inputdir"Sorted_bam/"${basenm}"_sorted.bam"
	
	# index the bam file
	samtools index $inputdir"Sorted_bam/"${basenm}"_sorted.bam"
done


for sf in $(ls $sorted_dir/*_sorted.bam); do
	bw_nam=$(echo ${sf} | sed 's/\_sorted.bam//' |cut -d "/" -f 2)
	echo "Coverage for " $bw_nam " in progress ..."
	bamCoverage -b $sf -o $outdir${bw_nam}.bw -of "bigwig" -p 5  --ignoreForNormalization chrX chrM -bs 20 --smoothLength 60  --normalizeUsingRPKM
done

