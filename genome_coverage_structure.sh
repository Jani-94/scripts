##first deactivat conda env
##module load python
##current script is only for chr 16,17 and 18

tinycov covplot --out test.png -p2 -m 600 -s 10000 -r 20000 -w 16,17,18 pig_merge.bam  

##Variant calling to test the allele ratio 
## -f think is right path to dir
## 
freebayes --bam Pig_IPECJ2_chr17.bam --use-best-n-alleles 4 -f ../../../../shared/public_data_store/genomes/pig/Sscrofa11.1.fa --min-base-quality 20 --min-mapping-quality 30 --min-alternate-fraction 0.2 --haplotype-length 0 --min-alternate-count 3 | vcffilter -f 'QUAL > 20' | bgzip -c > Pig_IPECJ2_chr17.vcf.gz

##Check the heterozygous SNPs and make histogram like this:
zcat Pig_IPECJ2.vcf.gz | cut -f10 | grep "0/1:" | cut -d':' -f2,3 | cut -d',' -f1 | sed 's/:/\t/g' | awk '{print $2/$1}' | awk -v size=0.025  '{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }'

zcat ../SNV_calling/snv_freebayes_2.vcf.gz | awk '$1 == "16"' | cut -f10 | grep "0/1:" | cut -d':' -f2,3 | cut -d',' -f1 | sed 's/:/\t/g' | awk '{print $2/$1}' | awk -v size=0.025  '{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax;bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }'

