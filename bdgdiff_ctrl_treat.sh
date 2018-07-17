#!/usr/bin/bash

## Author:  Maria Kondili

## Apply Differential Bound detection among 2 samples representing different time-points or different conditions/treatments to be compared.

workdir="/home/kondmar/ChIPseq/bam/Trimmed/Filtered/MACS2_PeakCalling/"

cd $workdir

outdir_narrow=$workdir'Narrow/'
outdir_broad=$workdir'Broad/'


# --d1, --d2 : Find the tags after filtering in control, mentioned in _peaks.xls - MACS2 callpeak Output
tags_wt0=10080529
tags_wt120=10080529

wt0=$workdir"bdg/WT_000m_R1/WT_000m_R1_treat_pileup.bdg"
wt120=$workdir"bdg/WT_120m_R1/WT_120m_R1_treat_pileup.bdg"

wt0_ctrl=$workdir"bdg/WT_000m_R1/WT_000m_R1_control_lambda.bdg"
wt120_ctrl=$workdir"bdg/WT_120m_R1/WT_120m_R1_control_lambda.bdg"


## Command in one-line :
 	
macs2 bdgdiff --t1 $wt0   --c1 $wt0_ctrl  --t2 $wt120  --c2 $wt120_ctrl  --d1 $tags_wt0 --d2  $tags_wt120  -l 120  --o-prefix  diffBind_0_vs_120


# l = min length of enriched region
 
## Outputs :
## Regions being enriched  in one of the conditions  : _cond1.bed  ,  _cond2.bed ,
## Regions having similar enrichment in both conditions :   _common.bed 



