##this one is right one!
bedtools bamtobed -i SL-29_ATAC_S13_markDupl.sorted_noMT.bam | awk -v OFS="\t" '{if ($6=="+") print $1,$2+4,$3+4; else if ($6=="-") print $1,$2-5,$3-5}' > SL-29_ATAC_S13_markDupl.sorted_noMT_shifted.bed
