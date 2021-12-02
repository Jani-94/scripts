#!/bin/bash

rsync -av /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/ChIP_seq/Chicken_crl/FINAL_January/peak_calling/H3K4me3_narrow_peaks.narrowPeak .
rsync -av /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/ChIP_seq/Chicken_crl/FINAL_January/peak_calling/H3K4me1_narrow_peaks.narrowPeak .
rsync -av /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/ChIP_seq/Chicken_crl/FINAL_January/peak_calling/H3K27ac_narrow_peaks.narrowPeak .
rsync -av /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/ChIP_seq/Chicken_crl/FINAL_January/peak_calling/H3K27me3_narrow_peaks.narrowPeak .
rsync -av /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/ChIP_seq/Chicken_crl/FINAL_January/peak_calling/CTCF_narrow_peaks.narrowPeak .

rsync -av /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/ChIP_seq/Pig_IPECJ2/REDO_analysis_FINAL_January/peak_calling/CTCF_narrow_peaks.narrowPeak Pig_CTCF_narrow_peaks.narrowPeak
rsync -av /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/ChIP_seq/Pig_IPECJ2/REDO_analysis_FINAL_January/peak_calling/H3K27me3_narrow_peaks.narrowPeak Pig_H3K27me3_narrow_peaks.narrowPeak
rsync -av /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/ChIP_seq/Pig_IPECJ2/REDO_analysis_FINAL_January/peak_calling/H3K27ac_narrow_peaks.narrowPeak Pig_H3K27ac_narrow_peaks.narrowPeak
rsync -av /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/ChIP_seq/Pig_IPECJ2/REDO_analysis_FINAL_January/peak_calling/H3K4me1_narrow_peaks.narrowPeak Pig_H3K4me1_narrow_peaks.narrowPeak
rsync -av /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/ChIP_seq/Pig_IPECJ2/REDO_analysis_FINAL_January/peak_calling/H3K4me3_narrow_peaks.narrowPeak Pig_H3K4me3_narrow_peaks.narrowPeak

### Run Homer chicken ###
#annotatePeaks.pl H3K4me3_narrow_peaks.narrowPeak Gallus_gallus.GRCg6a.dna.toplevel.fa -gtf Gallus_gallus.GRCg6a.103.chr.gtf > chicken_H3K4me3_output.tsv
#annotatePeaks.pl H3K4me1_narrow_peaks.narrowPeak Gallus_gallus.GRCg6a.dna.toplevel.fa -gtf Gallus_gallus.GRCg6a.103.chr.gtf > chicken_H3K4me1_output.tsv
#annotatePeaks.pl H3K27ac_narrow_peaks.narrowPeak Gallus_gallus.GRCg6a.dna.toplevel.fa -gtf Gallus_gallus.GRCg6a.103.chr.gtf > chicken_H3K27ac_output.tsv
#annotatePeaks.pl H3K27me3_narrow_peaks.narrowPeak Gallus_gallus.GRCg6a.dna.toplevel.fa -gtf Gallus_gallus.GRCg6a.103.chr.gtf > chicken_H3K27me3_output.tsv
#annotatePeaks.pl CTCF_narrow_peaks.narrowPeak Gallus_gallus.GRCg6a.dna.toplevel.fa -gtf Gallus_gallus.GRCg6a.103.chr.gtf > chicken_CTCF_output.tsv

### Run Homer pig ###
#annotatePeaks.pl Pig_H3K4me3_narrow_peaks.narrowPeak Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -gtf Sus_scrofa.Sscrofa11.1.103.chr.gtf > Pig_H3K4me3_output.tsv
#annotatePeaks.pl Pig_H3K4me1_narrow_peaks.narrowPeak Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -gtf Sus_scrofa.Sscrofa11.1.103.chr.gtf > Pig_H3K4me1_output.tsv
#annotatePeaks.pl Pig_H3K27ac_narrow_peaks.narrowPeak Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -gtf Sus_scrofa.Sscrofa11.1.103.chr.gtf > Pig_H3K27ac_output.tsv
#annotatePeaks.pl Pig_H3K27me3_narrow_peaks.narrowPeak Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -gtf Sus_scrofa.Sscrofa11.1.103.chr.gtf > Pig_H3K27me3_output.tsv
#annotatePeaks.pl Pig_CTCF_narrow_peaks.narrowPeak Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -gtf Sus_scrofa.Sscrofa11.1.103.chr.gtf > Pig_CTCF_output.tsv

### get methylation data and write to file ###

# Generate TSS file
#grep -P "\tgene\t" Gallus_gallus.GRCg6a.103.chr.gtf | awk '$7 == "+" {print $0}' | awk -F'[\t/="]' '{print $1"\t"$10"\t"$4-500"\t"$4+500}' > Chicken_TSS_500.bed
#grep -P "\tgene\t" Gallus_gallus.GRCg6a.103.chr.gtf | awk '$7 == "-" {print $0}' | awk -F'[\t/="]' '{print $1"\t"$10"\t"$5-500"\t"$5+500}' >> Chicken_TSS_500.bed

#grep -P "\tgene\t" Sus_scrofa.Sscrofa11.1.103.chr.gtf | awk '$7 == "+" {print $0}' | awk -F'[\t/="]' '{print $1"\t"$10"\t"$4-500"\t"$4+500}' > Pig_TSS_500.bed
#grep -P "\tgene\t" Sus_scrofa.Sscrofa11.1.103.chr.gtf | awk '$7 == "-" {print $0}' | awk -F'[\t/="]' '{print $1"\t"$10"\t"$5-500"\t"$5+500}' >> Pig_TSS_500.bed

#python methylation.py /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/int_analyses/Chicken_CRL/RNA_meth/rrbs_filtered_cg.bed.gz /lustre/nobackup/WUR/ABGC/vos198/analyses_FAANG_cell_lines/int_analyses/Chicken_CRL/RNA_meth/wgbs_filtered_cg.bed.gz Chicken_TSS_500.bed Gallus_gallus.GRCg6a.103.chr.gtf chicken
#python methylation.py /lustre/nobackup/WUR/ABGC/vos198/integrative_analysis/int_all_final/pig_rrbs_filtered_cg.bed.gz /lustre/nobackup/WUR/ABGC/vos198/integrative_analysis/int_all_final/pig_wgbs_filtered_cg.bed.gz Pig_TSS_500.bed Sus_scrofa.Sscrofa11.1.103.chr.gtf pig

### read chipseq data and combine with expression and methylation ###
#python chipseq_expression.py chicken_H3K4me3_output.tsv ../expression/rsem_expression_chicken.genes.results chicken
#python chipseq_expression.py chicken_H3K4me1_output.tsv ../expression/rsem_expression_chicken.genes.results chicken
#python chipseq_expression.py chicken_H3K27ac_output.tsv ../expression/rsem_expression_chicken.genes.results chicken
#python chipseq_expression.py chicken_H3K27me3_output.tsv ../expression/rsem_expression_chicken.genes.results chicken
#python chipseq_expression.py chicken_CTCF_output.tsv ../expression/rsem_expression_chicken.genes.results chicken

#python chipseq_expression.py Pig_H3K4me3_output.tsv ../expression/pig_calc_expr.genes.results pig
#python chipseq_expression.py Pig_H3K4me1_output.tsv ../expression/pig_calc_expr.genes.results pig
#python chipseq_expression.py Pig_H3K27ac_output.tsv ../expression/pig_calc_expr.genes.results pig
#python chipseq_expression.py Pig_H3K27me3_output.tsv ../expression/pig_calc_expr.genes.results pig
#python chipseq_expression.py Pig_CTCF_output.tsv ../expression/pig_calc_expr.genes.results pig

### concatenate output files ###
cat chicken_expr_chipseq_chicken_H3K4me3_output.tsv > chicken_chipseq_expression_methylation.tsv
cat chicken_expr_chipseq_chicken_H3K4me1_output.tsv chicken_expr_chipseq_chicken_H3K27ac_output.tsv chicken_expr_chipseq_chicken_H3K27me3_output.tsv | grep -v "PeakID" >> chicken_chipseq_expression_methylation.tsv
#chicken_expr_chipseq_chicken_CTCF_output.tsv 
#cat pig_expr_chipseq_Pig_H3K4me3_output.tsv > pig_chipseq_expression_methylation.tsv
#cat pig_expr_chipseq_Pig_H3K4me1_output.tsv pig_expr_chipseq_Pig_H3K27ac_output.tsv pig_expr_chipseq_Pig_H3K27me3_output.tsv pig_expr_chipseq_Pig_CTCF_output.tsv | grep -v "PeakID" >> pig_chipseq_expression_methylation.tsv

## generate plots
python generate_plots_geneswitch.py chicken_chipseq_expression_methylation.tsv chicken rsem_expression_chicken.genes.results 
#python generate_plots_geneswitch.py pig_chipseq_expression_methylation.tsv pig ../expression/pig_calc_expr.genes.results

## courgette paprika komkommer ui champignon