#!/bin/bash

## Author : Maria Kondili 
## Subject : Find known motifs over the enriched regions 
## Further Info on the algorithm : Doc: http://homer.ucsd.edu/homer/ngs/peakMotifs.html  
cd /home/kondmar/Documents/NV/IRF5_ChIPseq/MACS2_PeakCalling/

output_dir="/home/kondmar/Documents/NV/IRF5_ChIPseq/HOMER_MotifDiscovery/"
peaks_WT1="WT30m_broadPeak.bed"
peaks_WT2="WT120m_broadPeak.bed"

	
## De novo motif Discovery with given length 

##  Motif analysis with Mask and background size
findMotifsGenome.pl $peaks_WT1 mm10 $output_dir'WT_30m/' -len 12 -size 100 -mask
findMotifsGenome.pl $peaks_WT2 mm10 $output_dir'WT_120m/' -len 12 -size 100 -mask




