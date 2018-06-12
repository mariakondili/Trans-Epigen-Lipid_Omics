#!/bin/bash
mkdir /media/kondmar/my2TB/PROJECTS/NV/IRF5_ChIPseq/Bigwigs/
cd /media/kondmar/my2TB/PROJECTS/NV/IRF5_ChIPseq/bam/Bowtie2_Trimmed/WT/

#for f in $(ls chip/*.bam); do

#	echo $f
#	basenm=$(echo ${f} | sed 's/.bam//')
	# Sort BAM file
#	samtools sort -O BAM $f -o "Sorted_bam/"${basenm}"_sorted.bam"
	# index the bam file
#	samtools index ${basenm}_sorted.bam
#done


for sf in $(ls *_sorted.bam); do
	bw_nam=$(echo ${sf} | sed 's/\_sorted.bam//')
	echo "Coverage for " $bw_nam " will be calculated now.."
	bamCoverage -b $sf -o ../../../Bigwigs/${bw_nam}.bw -of "bigwig" -p 5  --ignoreForNormalization chrX chrM -bs 20 --smoothLength 60  --normalizeUsingRPKM
	
done


##### In EMTAB-2661 
##> Use the bam files from IRF5_ChIPseq of Mano's directory

source_dir='/media/Timecapsule/Presentations/Mano/IRF5/IRF5ChIP-Seq/IRF5_ChIP_BMDM/BAM'
mkdir Sorted_bam/
for f in $(ls $source_dir/*.bam); do
	
	basenm=`echo $f |rev | cut -d"/" -f 1 |rev | sed 's/.bam//'`  ## rev = reads backwards, from end to beginning, so 1st field is the file.bam
	# Sort BAM file
	echo $basenm
	#samtools sort -O BAM $f -o "Sorted_bam/"${basenm}"_sorted.bam"
	# index the bam file
	samtools index "Sorted_bam/"${basenm}"_sorted.bam"
done


for sf in $(ls Sorted_bam/*_sorted.bam); do
	bw_nam=$(echo ${sf} | sed 's/\_sorted.bam//' |cut -d "/" -f 2)
	echo "Coverage for " $bw_nam " will be calculated now.."
	bamCoverage -b $sf -o Bigwigs/${bw_nam}.bw -of "bigwig" -p 5  --ignoreForNormalization chrX chrM -bs 20 --smoothLength 60  --normalizeUsingRPKM
done


cd /media/kondmar/my2TB/PROJECTS/NV/IRF5_ChIPseq/bam/BAM_Trim_q10_RmDup/KO/ 
bamCoverage -b KO_000m_R1_q10_rmDup.bam -o ../../../Bigwigs/KO_000m_R1.bw -of "bigwig" -p 5  --ignoreForNormalization chrX chrM -bs 20 --smoothLength 60  --normalizeUsingRPKM


