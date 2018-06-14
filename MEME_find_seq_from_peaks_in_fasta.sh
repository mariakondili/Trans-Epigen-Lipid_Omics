#!/bin/bash

## Author : Maria Kondili

###Find DNA sequence of peaks --> to be sent to MEME for motif analysis

# fasta Ref to use : 
ref_genome='/home/kondmar/GenomeRef/GRCm38.p5.genome.fa' # Should have chrom as header per sequence in the fasta ( mm10.ref.fa is not structured with header)
#Index 
samtools faidx $ref_genome  ## Creates Output .fa.fai in same Directoyry as Input

#> bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf> -fo <fasta> 

##> Bed files to be used : 
bed_dir='/home/kondmar/ChIPseq/MACS2_PeakCalling/Broad/'
bed_0m=$bed_dir'MergedBAMs/WT0m_mrg_broadPeaks.bed'
bed_120m=$bed_dir'MergedBAMs/WT120m_mrg_broadPeaks.bed'
bed_30m=$bed_dir'WT30m/WT30m_broadPeaks.bed'

out_dir=$bed_dir'MEME_MOTIFS/'

bedtools getfasta -fi $ref_genome -bed $bed_0m -fo $out_dir'seq_in_pks_0m.fa' 
bedtools getfasta -fi $ref_genome -bed $bed_120m -fo $out_dir'seq_in_pks_120m.fa' 
bedtools getfasta -fi $ref_genome -bed $bed_30m -fo $out_dir'seq_in_pks_30m.fa'

##WARNING:  chromosome (chrUn_JH584304) was not found in the FASTA file

