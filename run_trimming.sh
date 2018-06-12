#!/bin/bash

## Author : Maria Kondili 

## Subject : Run Trimmomatic for trimming adapters from fastq files

# java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz 
# 						   	 output_forward_paired.fq.gz output_forward_unpaired.fq.gz 
# 						   	 output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz 
# 						   	 ILLUMINACLIP:/root/path/TruSeq3-PE-2.fa:2:30:10 
# 						   	 LEADING:3 TRAILING:3 
# 						   	 SLIDINGWINDOW:5:20
# 						   	 MINLEN:50

#> ILLUMINACLIP:TruSeq3-PE-2.fa: A:M:K 
# TruSeq3-PE.fa : for sequencing done in HiSeq. Should be downloaded and the full path given along with fasta filename
# A=Maximum mismatch count which will still allow a full match to be performed
# M=How accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment
# K=How accurate the match between any adapter etc. sequence must be against a read

trim_path="/home/Tools/Trimmomatic/Trimmomatic-0.36/" # to be installed from 
input_dir="/home/Projects/ChIPseq/fastq/"
output_dir=$input_dir'Trimming/Fastq_trimmed/'


##> PE: Paired-end sequenced fastq files 
##  Suffix : _R1-1.fastq.gz   _R1-2.fastq.gz
java -jar $trim_path'trimmomatic-0.36.jar' PE -threads 8 -phred33  $input_dir"Sample_R1-1.fastq.gz" $input_dir"Sample_R1-2.fastq.gz" $output_dir"Sample_R1_fwd_paired.fq.gz"   WT-0m_R1_fwd_unpaired.fq.gz $output_dir"Sample_R1_rvs_paired.fq.gz"  WT-0m_R1_rvs_unpaired.fq.gz  ILLUMINACLIP:$trim_path/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3  SLIDINGWINDOW:4:20 MINLEN:50 


##> SE:  Single-end sequenced fastq files
outdir_30="/media/kondmar/my2TB/PROJECTS/NV/IRF5_ChIPseq/fastq/Trimming/Fastq_trimmed/"
java -jar $trim_path'trimmomatic-0.36.jar' SE -threads 8  -phred33 $input_dir"SE_Sample_R1.fastq.gz"  $outdir_30"SE_Sample_R1_trimmed.fq.gz" ILLUMINACLIP:$trim_path/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3  SLIDINGWINDOW:4:20  MINLEN:50
