import seaborn as sns
import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr  
import numpy as np
import subprocess 

df = pd.read_csv(sys.argv[1], sep='\t') ## pig
#print(pearsonr(df['Peak Score'],df['TPM']))
df['TPM'] = np.log2(df['TPM']) 
df['Peak Score'] = np.log2(df['Peak Score']) 

## Final heatmap ##
df_new = df.filter(["Mark","TPM","Peak Score","RRBS_TSS", "WGBS_TSS","RRBS_GB","WGBS_GB"], axis=1)
dataCorr = df_new.groupby('Mark')[["TPM","Peak Score","RRBS_TSS", "WGBS_TSS","RRBS_GB","WGBS_GB"]].corr()#.iloc[0::2][['Peak Score']]

dataCorr = dataCorr.drop('Peak Score', axis=1)
dataCorr = dataCorr.drop(labels='RRBS_GB', level=1)
dataCorr = dataCorr.drop(labels='RRBS_TSS', level=1)
dataCorr = dataCorr.drop(labels='WGBS_GB', level=1)
dataCorr = dataCorr.drop(labels='WGBS_TSS', level=1)
dataCorr = dataCorr.drop(labels='TPM', level=1)

## Within group correlation ##
# df_new = df.filter(["Class","Mark","TPM","Peak Score","RRBS_TSS", "WGBS_TSS","RRBS_GB","WGBS_GB"], axis=1)
# dataCorr = df_new.groupby('Class')[["TPM","Peak Score","RRBS_TSS", "WGBS_TSS","RRBS_GB","WGBS_GB"]].corr()#.iloc[0::2][['Peak Score']]
# print(dataCorr)

# dataCorr = dataCorr.drop('TPM', axis=1)
# dataCorr = dataCorr.drop(labels='RRBS_GB', level=1)
# dataCorr = dataCorr.drop(labels='RRBS_TSS', level=1)
# dataCorr = dataCorr.drop(labels='WGBS_GB', level=1)
# dataCorr = dataCorr.drop(labels='WGBS_TSS', level=1)
# dataCorr = dataCorr.drop(labels='Peak Score', level=1)

## Expression vs methylation and peak score
#dataCorr = df.filter(["TPM","Peak Score","RRBS_TSS", "WGBS_TSS","RRBS_GB","WGBS_GB"], axis=1).corr()
#print(dataCorr)

#dataCorr = dataCorr.drop(labels='RRBS_GB', level=1))

# groups = ['TPM','Peak Score','RRBS_TSS', 'WGBS_TSS','RRBS_GB','WGBS_GB']
# df2 = pd.DataFrame()
# for i in range( len(groups)-1): 
    # df2 = df2.append( df_new.groupby('Mark')[groups].corr().stack()
                        # .loc[:,groups[i],groups[i+1]:].reset_index() )

# df2.columns = ['ID', 'v1', 'v2', 'corr']
# df2.set_index(['ID','v1','v2']).sort_index()
# print(df2)
# df_m = df2.groupby(["ID", "v1", 'corr']).size().unstack(level=0)
#dataCorr = dataCorr.mask(np.tril(np.ones(dataCorr.shape)).astype(np.bool))
#print(dataCorr)


#corr = df_new.groupby('Mark')[["TPM","Peak Score","RRBS_TSS", "WGBS_TSS","RRBS_GB","WGBS_GB"]].corr()#.iloc[0::2][['Peak Score']]
#print(corr)

#plt.figure(0)
#sns_plot_chicken = sns.scatterplot(x="Peak Score", y="TPM", log_scale=True, data=df) #, showfliers = False, order=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','30','31','32','33','Z'])
#splot = sns.regplot(x="Peak Score", y="TPM", 
#                    data=df, scatter_kws={'alpha':0.15}, fit_reg=True, scatter=True)
#splot.set(xscale="log")
#splot.set(yscale="log")
# plt.xticks(rotation=90)
#splot.figure.savefig(sys.argv[1]+"_chicken_scatter_expr_chip.png")

#plt.figure(1)
#ax = sns.violinplot([np.log2(vlow),np.log2(low),np.log2(medium),np.log2(high),np.log2(vhigh)])
#ax.savefig(sys.argv[1]+"_chicken_violin_expr_chip.png")
# plt.clf()

#plt.figure(2)
plt.figure(0)
fig_dims = (20, 8)
fig, ax = plt.subplots(figsize=fig_dims)
splot = sns.boxplot(x="Class", y="Peak Score", hue="Mark", data=df, order=['TPM<1','TPM<5','TPM<20','TPM<100','TPM>=100'], hue_order=["H3K4me3","H3K27ac","H3K4me1","H3K27me3"])
#splot.set(ylim=(4, 13))
splot.set_xlabel("Expression Class",fontsize=15,labelpad=10)
splot.set_ylabel("Peak Score",fontsize=15)
splot.tick_params(labelsize=12)
splot.legend(fontsize='12', title_fontsize='12')
plt.savefig(sys.argv[2]+"_boxplot_box_expr_chip.png")

plt.figure(1)
fig_dims = (20, 8)
fig, ax = plt.subplots(figsize=fig_dims)
splot = sns.violinplot(x="Class", y="Peak Score", hue="Mark", data=df, order=['TPM<1','TPM<5','TPM<20','TPM<100','TPM>=100'], hue_order=["H3K4me3","H3K27ac","H3K4me1","H3K27me3"])
splot.set_xlabel("Expression Class",fontsize=15,labelpad=10)
splot.set_ylabel("Peak Score",fontsize=15)
splot.tick_params(labelsize=12)
splot.legend(fontsize='12', title_fontsize='12')
plt.savefig(sys.argv[2]+"_violin_box_expr_chip.png")

## Methylation ##
df_plot = df.melt(id_vars='Class', var_name='dataset', value_vars=["RRBS_TSS", "WGBS_TSS","RRBS_GB","WGBS_GB"])

plt.figure(2)
fig_dims = (20, 8)
fig, ax = plt.subplots(figsize=fig_dims)
splot = sns.boxplot(x="Class", y="value", hue='dataset', data=df_plot, order=['TPM<1','TPM<5','TPM<20','TPM<100','TPM>=100'])
splot.set_xlabel("Expression Class",fontsize=15,labelpad=10)
splot.set_ylabel("Methylation level",fontsize=15)
splot.tick_params(labelsize=12)
splot.legend(loc='upper right',fontsize='12', title_fontsize='12')
plt.savefig(sys.argv[2]+"_box_expr_meth.png")

plt.figure(3)
fig_dims = (20, 8)
fig, ax = plt.subplots(figsize=fig_dims)
splot = sns.violinplot(x="Class", y="value", hue='dataset', data=df_plot, order=['TPM<1','TPM<5','TPM<20','TPM<100','TPM>=100'])
splot.set_xlabel("Expression Class",fontsize=15,labelpad=10)
splot.set_ylabel("Methylation level",fontsize=15)
splot.tick_params(labelsize=12)
splot.legend(loc='upper right', fontsize='12', title_fontsize='12')
plt.savefig(sys.argv[2]+"_violin_expr_meth.png")

plt.figure(4)
fig_dims = (20, 8)
fig, ax = plt.subplots(figsize=fig_dims)
splot = sns.heatmap(dataCorr, annot=True, annot_kws={"fontsize":12})
splot.set_ylabel("")
splot.tick_params(labelsize=12)
#splot.set(ylim=(4, 13))
plt.savefig(sys.argv[2]+"_heat_expr_chip.png")

# sns_plot_pig = sns.boxplot(x="Chr", y="TPM", data=df, showfliers = False, order=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','X'])
# sns_plot_pig.figure.savefig("pig_boxplot_expr_per_chr.png")

## test ## 
# Read methylation data and fill dictionary  
outfile_expression_methylation = open(sys.argv[2]+"_expr_meth.tsv","w")
outfile_expression_methylation.write("TPM\tRRBS_TTS\tWGBS_TSS\tRRBS_GB\tWGBS_GB\n")
meth_dic={}
chicken_meth = open(sys.argv[2]+"_methylation.tsv","r")
for meth in chicken_meth:
    meth=meth.split("\t")
    meth_dic[meth[0]] = [meth[1],meth[2],meth[3],meth[4].strip()]

expr_dic={}
chicken_expr = open(sys.argv[3],"r")
header = chicken_expr.readline()
for expr in chicken_expr:
    if expr.split("\t")[0] in meth_dic:
        expr_dic[expr.split("\t")[0]] = expr.split("\t")[5] ## dictionary with gene symbol (key) and TPM expression
        outlist=[expr.split("\t")[5]]+meth_dic[expr.split("\t")[0]] 
        outfile_expression_methylation.write("\t".join(outlist)+"\n")
outfile_expression_methylation.close()    

df1 = pd.read_csv(sys.argv[2]+"_expr_meth.tsv", sep='\t')
df1['TPM'] = np.log2(df1['TPM'])
dataCorr=df1.corr()
plt.figure(5)
fig_dims = (20, 8)
fig, ax = plt.subplots(figsize=fig_dims)
splot = sns.heatmap(dataCorr, annot=True, annot_kws={"fontsize":12})
splot.tick_params(labelsize=12)
#splot.set(ylim=(4, 13))
plt.savefig(sys.argv[2]+"_heat_expr_meth.png")
#print(pearsonr(df['Peak Score'],df['TPM']))
#df['TPM'] = np.log2(df['TPM']) 
#df['Peak Score'] = np.log2(df['Peak Score']) 
