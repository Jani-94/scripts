
import seaborn as sns
import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr  
import numpy as np
import subprocess 

## combine expression and chipseq
def combine_expression_methylation_chipseq():

    # Read expression data and fill dictionary
    expr_dic={}
    chicken_expr = open(sys.argv[2],"r")
    header = chicken_expr.readline()
    for expr in chicken_expr:
        expr_dic[expr.split("\t")[0]] = expr.split("\t")[5] ## dictionary with gene symbol (key) and TPM expression

    # Read methylation data and fill dictionary  
    meth_dic={}
    chicken_meth = open(sys.argv[3]+"_methylation.tsv","r")
    for meth in chicken_meth:
        meth=meth.split("\t")
        meth_dic[meth[0]] = [meth[1],meth[2],meth[3],meth[4]]

    # Read chipseq data and output file with chipseq info, expression, and methylation
    vlow,low,medium,high,vhigh=[],[],[],[],[]
    chicken_chipseq = open(sys.argv[1],"r")
    header = chicken_chipseq.readline()
    outfile=open(sys.argv[3]+"_expr_chipseq_"+sys.argv[1],"w")
    outfile.write("PeakID (cmd=annotatePeaks.pl H3K4me3_narrow_peaks.narrowPeak Gallus_gallus.GRCg6a.dna.toplevel.fa -gtf Gallus_gallus.GRCg6a.103.chr.gtf)	Chr	Start	End	Strand	Peak Score	Focus Ratio/Region Size	Annotation	Detailed Annotation	Distance to TSS	Nearest PromoterID	Entrez ID	Nearest Unigene	Nearest Refseq	Nearest Ensembl	Gene Name	Gene Alias	Gene Description	Gene Type\tTPM\tSignal\tClass\tMark\tRRBS_TSS\tWGBS_TSS\tRRBS_GB\tWGBS_GB\n")

    #list_features=["promoter-TSS", "intron"]
    list_features=["promoter-TSS"]
    for chip in chicken_chipseq:
        chip_split=chip.split("\t")
        if any(s in chip_split[8] for s in list_features):
            #print(chip_split)
            geneid=chip_split[11]
            expression = float(expr_dic[geneid])
            rrbs_tss,wgbs_tss,rrbs_gb,wgbs_gb = meth_dic[geneid]
            peak_id=chip_split[0]
            signal_value = "0"
            classvalue=""
            #signal_value = subprocess.getoutput("grep -P '"+peak_id+"\t' H3K27ac_narrow_peaks.narrowPeak | cut -f10")
            if expression < 1.0:
                classvalue = "TPM<1"
            elif expression < 5.0:
                classvalue = "TPM<5"
            elif expression < 20.0:
                classvalue = "TPM<20"
            elif expression < 100:
                classvalue = "TPM<100"
            else:
                classvalue = "TPM>=100"

            outfile.write(chip.strip()+"\t"+str(expression)+"\t"+signal_value+"\t"+classvalue+"\t"+sys.argv[1].split("_")[1]+"\t"+rrbs_tss+"\t"+wgbs_tss+"\t"+rrbs_gb+"\t"+wgbs_gb.strip()+"\n") ## write homer results + expression value
    outfile.close()    
        
    # df = pd.read_csv("chicken_expr_chipseq_"+sys.argv[1], sep='\t') ## pig
    # #print(pearsonr(df['Peak Score'],df['FPKM']))
    # df['FPKM'] = np.log2(df['FPKM']) 
    # df['Peak Score'] = np.log2(df['Peak Score']) 
    # df_new = df.filter(['Mark','FPKM','Peak Score'], axis=1)
    # corr = df_new.groupby('Mark').corr()

    #plt.figure(0)
    #sns_plot_chicken = sns.scatterplot(x="Peak Score", y="TPM", log_scale=True, data=df) #, showfliers = False, order=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','30','31','32','33','Z'])
    #splot = sns.regplot(x="Peak Score", y="FPKM", 
    #                    data=df, scatter_kws={'alpha':0.15}, fit_reg=True, scatter=True)
    #splot.set(xscale="log")
    #splot.set(yscale="log")
    # plt.xticks(rotation=90)
    #splot.figure.savefig(sys.argv[1]+"_chicken_scatter_expr_chip.png")

    #plt.figure(1)
    #ax = sns.violinplot([np.log2(vlow),np.log2(low),np.log2(medium),np.log2(high),np.log2(vhigh)])
    #ax.savefig(sys.argv[1]+"_chicken_violin_expr_chip.png")
    # plt.clf()

    # #plt.figure(2)
    # plt.figure(0)
    # fig_dims = (18, 8)
    # fig, ax = plt.subplots(figsize=fig_dims)
    # splot = sns.boxplot(x="Class", y="RRBS_TSS", hue="Mark", data=df, order=['vlow','low','medium','high','vhigh'], hue_order=["H3K4me3","H3K27ac","H3K4me1","H3K27me3","CTCF"])
    # #splot.set(ylim=(4, 13))
    # plt.savefig(sys.argv[1]+"_chicken_box_expr_chip.png")

    # plt.figure(1)
    # fig_dims = (18, 8)
    # fig, ax = plt.subplots(figsize=fig_dims)
    # splot = sns.heatmap(corr)
    # #splot.set(ylim=(4, 13))
    # plt.savefig(sys.argv[1]+"_chicken_heat_expr_chip.png")

    # sns_plot_pig = sns.boxplot(x="Chr", y="TPM", data=df, showfliers = False, order=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','X'])
    # sns_plot_pig.figure.savefig("pig_boxplot_expr_per_chr.png")

#get_and_write_methylation()
combine_expression_methylation_chipseq()