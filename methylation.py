import sys
import subprocess

def get_m_level(chr,start,stop): ## include strand
    rrbs=subprocess.getoutput("/cm/shared/apps/htslib/gcc/64/1.9/bin/tabix "+sys.argv[1]+" "+str(chr)+":"+str(start)+"-"+str(stop)+" | awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }'")
    wgbs=subprocess.getoutput("/cm/shared/apps/htslib/gcc/64/1.9/bin/tabix "+sys.argv[2]+" "+str(chr)+":"+str(start)+"-"+str(stop)+" | awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }'")  
    return(rrbs, wgbs)

def get_and_write_methylation():    
    TSS_coordinates = open(sys.argv[3], "r")
    GTF_annotation = open(sys.argv[4],"r")

    dic_tss={} ## create empty dic
    for gene_c in TSS_coordinates: ## Loop file
        gene_c=gene_c.strip().split()
        dic_tss[gene_c[1]] = [gene_c[0],gene_c[2],gene_c[3]] ## put strand colum in value list

    dic_gb={} ## create empty dic
    for gene in GTF_annotation: ## Loop file
        if gene.startswith("#"):
            continue
        gene=gene.split("\t")
        if gene[2] != "gene":
            continue
        gene_id=gene[8].split(" ")[1].strip(';').strip('"')
        chr,start,end = gene[0],gene[3],gene[4]
        dic_gb[gene_id] = [chr,start,end] ## put strand colum in value list

    meth_output=open(sys.argv[5]+"_methylation.tsv","w")
    for gene in dic_gb:
        chr, start, stop = dic_gb[gene] ## get GB coordinates
        gb_mlevel_rrbs,gb_mlevel_wgbs = get_m_level(chr,int(start),int(stop))
        chr, start, stop = dic_tss[gene] ## get Tss coordinates
        tss_mlevel_rrbs,tss_mlevel_wgbs = get_m_level(chr,int(start),int(stop)) 
        table = [gene,tss_mlevel_rrbs,tss_mlevel_wgbs,gb_mlevel_rrbs,gb_mlevel_wgbs]
        meth_output.write("\t".join(table)+"\n")   

get_and_write_methylation()